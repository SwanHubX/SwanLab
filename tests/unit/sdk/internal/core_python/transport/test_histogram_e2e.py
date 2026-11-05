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

import swanlab
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import RunState
from swanlab.sdk.cmd import sync as sync_cmd
from swanlab.sdk.internal.core_python import client as core_client
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.utils import PrepareExperimentStartResult
from swanlab.sdk.internal.settings import Settings


class _PlaceholderBackend:
    """本地 HTTP 测试后端，接收指标和列定义，并提供实验摘要。"""

    def __init__(self):
        self.metrics_requests: list = []
        self.column_requests: list = []
        self.summary: dict = {"log": None, "media": None, "scalar": None, "histogram": None}
        self.summary_requests = 0
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
                elif self.path.endswith("/series"):
                    backend.column_requests.append(data)
                    self._send("SUCCESS")
                else:
                    self._send({"code": 404, "message": "not found"}, status=404)

            def do_GET(self):
                if "/house/metrics/summaries/" in self.path:
                    backend.summary_requests += 1
                    self._send(backend.summary)
                else:
                    self._send({"code": 404, "message": "not found"}, status=404)

            def do_PUT(self):
                length = int(self.headers.get("Content-Length", 0))
                self.rfile.read(length)
                if self.path.endswith("/profile"):
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
    assert metric["data"] == {"values": [2, 5, 0], "bins": [-1.0, 0.0, 0.5, 2.0]}
    assert "create_time" in metric


def test_upload_histogram_empty_batch_skips_http(placeholder_client, tmp_path: Path):
    """空批次不发出任何请求。"""
    sender = _make_sender(tmp_path)
    sender.upload("histogram", [])
    assert placeholder_client.metrics_requests == []


def test_upload_histogram_splits_http_requests(placeholder_client, tmp_path: Path):
    sender = _make_sender(tmp_path)
    records = [_make_histogram_record(step=step) for step in range(513)]

    result = sender.upload("histogram", records)

    requests = placeholder_client.metrics_requests
    assert [len(request["metrics"]) for request in requests] == [512, 1]
    assert all(request["type"] == "histogram" for request in requests)
    metrics = [metric for request in requests for metric in request["metrics"]]
    assert [metric["index"] for metric in metrics] == list(range(513))
    assert all("epoch" not in metric for metric in metrics)
    assert result.confirmed_records == records
    assert result.remaining_records == []
    assert result.failures == []


@pytest.mark.parametrize(
    "new_experiment, summary_step, range_start, expected_steps",
    [
        (True, None, None, [0, 1, 2]),
        (False, None, None, [0, 1, 2]),
        (False, 1, None, [2]),
        (False, 2, None, []),
        (False, 1, 0, [0, 1, 2]),
        (False, 1, 1, [1, 2]),
    ],
)
def test_sync_offline_histogram_file_to_http(
    placeholder_client, tmp_path: Path, monkeypatch, new_experiment, summary_step, range_start, expected_steps
):
    run = swanlab.init(mode="offline", project="histogram-sync", id="histogram-sync", log_dir=str(tmp_path))
    run_dir = run._ctx.run_dir
    for step in range(3):
        swanlab.log({"weights": swanlab.Histogram(np_histogram=([2, 1], [0.0, 1.0, 2.0]))}, step=step)
    swanlab.finish()

    reader = DataStoreReader()
    reader.open(next(run_dir.glob("run-*.swanlab")))
    try:
        source = []
        for data in reader:
            record = Record.FromString(data)
            if record.HasField("histogram"):
                source.append(record.histogram)
    finally:
        reader.close()
    assert [record.step for record in source] == [0, 1, 2]
    if summary_step is not None:
        placeholder_client.summary["histogram"] = [{"key": "weights", "step": summary_step}]

    prepared = PrepareExperimentStartResult(
        username="alice",
        project="demo",
        name="histogram-sync",
        color="#abcdef",
        project_data={
            "cuid": "project-id",
            "name": "demo",
            "version": 1,
            "group": {"username": "alice"},
            "path": "/alice/demo",
            "visibility": "PRIVATE",
            "_count": {"experiments": 0, "contributors": 0, "collaborators": 0, "clones": 0},
        },
        experiment={
            "cuid": "experiment-id",
            "slug": "histogram-sync",
            "name": "histogram-sync",
            "createdAt": "2024-08-01T00:00:00Z",
        },
        new_experiment=new_experiment,
    )
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.prepare_experiment_start", lambda record: prepared)
    stopped = []
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.stop_experiment", lambda *args, **kwargs: stopped.append(kwargs["state"])
    )

    sync_cmd.sync(
        run_dir,
        settings=Settings(run=Settings.Run(id="histogram-sync")),
        start_step=range_start,
    )

    requests = [request for request in placeholder_client.metrics_requests if request["type"] == "histogram"]
    metrics = [metric for request in requests for metric in request["metrics"]]
    assert [metric["index"] for metric in metrics] == expected_steps
    assert all(metric["key"] == "weights" for metric in metrics)
    assert all(metric["data"] == {"values": [2, 1], "bins": [0.0, 1.0, 2.0]} for metric in metrics)
    assert all("epoch" not in metric for metric in metrics)
    assert placeholder_client.summary_requests == (0 if new_experiment else 1)
    if summary_step is None:
        columns = [column for request in placeholder_client.column_requests for column in request["series"]]
        assert any(column["key"] == "weights" and column["type"] == "HISTOGRAM" for column in columns)
    assert stopped == [RunState.RUN_STATE_FINISHED]
