"""
@author: mx
@file: test_histogram_e2e.py
@time: 2026/10/28
@description: 直方图 HTTP 集成测试：sender 编码 → client HTTP → 本地 HTTP 服务
"""

import json
import threading
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python import client as core_client
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender


class _PlaceholderBackend:
    """占位后端：实现 login 与 /house/metrics 两个端点，记录收到的指标请求体"""

    def __init__(self):
        self.metrics_requests: list = []
        backend = self

        class Handler(BaseHTTPRequestHandler):
            def _send(self, payload, status=200):
                body = json.dumps(payload).encode("utf-8")
                self.send_response(status)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_POST(self):
                length = int(self.headers.get("Content-Length", 0))
                data = json.loads(self.rfile.read(length) or b"{}")
                if self.path.endswith("/login/api_key"):
                    expired = (datetime.now(timezone.utc) + timedelta(days=30)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                    self._send({"sid": "test-sid", "expiredAt": expired, "userInfo": {"username": "tester"}})
                elif self.path.endswith("/house/metrics"):
                    backend.metrics_requests.append(data)
                    self._send("SUCCESS")
                else:
                    self._send({"code": 404, "message": "not found"}, status=404)

            def log_message(self, format: str, *args: object) -> None:
                pass

        self._server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)

    def start(self) -> str:
        self._thread.start()
        host, port = self._server.server_address[:2]
        return f"http://{host}:{port}"

    def stop(self) -> None:
        self._server.shutdown()
        self._server.server_close()
        self._thread.join(timeout=5)


@pytest.fixture
def placeholder_client():
    """启动占位后端并把全局运行时客户端指向它，测试结束后还原全局状态。"""
    backend = _PlaceholderBackend()
    base_url = backend.start()
    leaked = core_client.exists()
    if leaked:
        core_client.reset()
    core_client.new(api_key="test-key", base_url=base_url)
    try:
        yield backend
    finally:
        core_client.reset()
        backend.stop()


def _make_sender(tmp_path: Path) -> HttpRecordSender:
    ctx = CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=tmp_path,
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
        )
    )
    ctx.set_online_params(
        username="alice",
        project="demo",
        project_id="project-id",
        project_version=1,
        experiment_id="experiment-id",
    )
    return HttpRecordSender(ctx=ctx)


def _make_histogram_record(key: str = "train/w", step: int = 3) -> Record:
    ts = Timestamp()
    ts.GetCurrentTime()
    return Record(
        histogram=HistogramRecord(
            key=key,
            step=step,
            type=ColumnType.COLUMN_TYPE_HISTOGRAM,
            timestamp=ts,
            value={"values": [2, 5, 0], "bins": [-1.0, 0.0, 0.5, 2.0]},
        )
    )


def test_upload_histogram_http_roundtrip(placeholder_client, tmp_path: Path):
    """sender.upload 分发 histogram 批次 → 占位后端收到 House 写入协议规定的请求体。"""
    sender = _make_sender(tmp_path)
    sender.upload("histogram", [_make_histogram_record(), _make_histogram_record(key="train/b", step=4)])

    assert len(placeholder_client.metrics_requests) == 1
    payload = placeholder_client.metrics_requests[0]
    assert payload["type"] == "histogram"
    assert payload["projectId"] == "project-id"
    assert payload["experimentId"] == "experiment-id"
    assert len(payload["metrics"]) == 2
    assert all("epoch" not in metric for metric in payload["metrics"])
    metric = payload["metrics"][0]
    assert metric["key"] == "train/w"
    assert metric["index"] == 3
    assert metric["data"] == {"_type": "histogram", "values": [2, 5, 0], "bins": [-1.0, 0.0, 0.5, 2.0]}
    assert "create_time" in metric


def test_upload_histogram_empty_batch_skips_http(placeholder_client, tmp_path: Path):
    """空批次不发出任何请求。"""
    sender = _make_sender(tmp_path)
    sender.upload("histogram", [])
    assert placeholder_client.metrics_requests == []
