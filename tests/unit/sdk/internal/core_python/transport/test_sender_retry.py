"""
重试设计契约测试：验证上传链路中 ApiError 的传播规则和 retries=0 约定。

重试设计分三层：
  1. API 层（api/upload.py）：所有 /house/metrics 和 column 上传使用 retries=0，
     HTTP 层不重试，重试交由 transport 层处理。
  2. sender.upload（sender.py）：重试决定点。5xx ApiError → re-raise（transport 重试）；
     4xx ApiError → swallow（业务拒绝，不重试）。
  3. Dispatch（dispatch.py）：safe.block 包裹 sender.upload，捕获所有异常返回 (False, records)。

契约：handler 不得将 HTTP API 调用包裹在 safe.block 中，
否则 sender.upload 永远看不到 ApiError，transport 层无法重试。
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord, ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.api import upload as upload_api
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker

# ── Helpers ──────────────────────────────────────────────────


class _FakeApiErrorResponse:
    """构造 ApiError 所需的最小 Response 替身。"""

    def __init__(self, status_code: int) -> None:
        self.status_code = status_code
        self.request = MagicMock()
        self.url = "https://api.test/upload"


def _make_api_error(status_code: int) -> ApiError:
    return ApiError(
        _FakeApiErrorResponse(status_code),
        method="POST",
        trace_id="trace-id",
        code="test-error",
        message=f"http {status_code}",
    )


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


def _make_scalar_record() -> Record:
    ts = Timestamp()
    ts.GetCurrentTime()
    return Record(
        scalar=ScalarRecord(
            key="train/loss",
            step=1,
            timestamp=ts,
            type=ColumnType.COLUMN_TYPE_SCALAR,
            value=ScalarValue(number=0.125),
        )
    )


def _make_log_record() -> Record:
    ts = Timestamp()
    ts.GetCurrentTime()
    return Record(
        log=LogRecord(
            line="training step done",
            level=LogLevel.LOG_LEVEL_INFO,
            timestamp=ts,
            epoch=1,
        )
    )


def _make_column_record() -> Record:
    return Record(column=ColumnRecord(column_key="metric-0", column_type=ColumnType.COLUMN_TYPE_SCALAR))


@pytest.fixture(autouse=True)
def _silence_console():
    """抑制 safe.block 的 console.trace 和 sender.upload 的 console.warning 输出。"""
    with (
        patch("swanlab.sdk.internal.pkg.console.trace"),
        patch("swanlab.sdk.internal.pkg.console.warning"),
    ):
        yield


# ── Section 1: retries=0 约定 ──────────────────────────────
# 所有上传 API 函数必须传 retries=0，HTTP 层不重试，重试交由 transport 层处理。


_RETRIES_ZERO_CASES = [
    pytest.param(
        upload_api.upload_scalar,
        ("pid", "eid"),
        {"metrics": [{"key": "loss", "index": 0, "data": 0.1, "create_time": "2024-01-01T00:00:00+08:00"}]},
        id="scalar",
    ),
    pytest.param(
        upload_api.upload_log,
        ("pid", "eid"),
        {"metrics": [{"level": "INFO", "epoch": 0, "message": "hi", "create_time": "2024-01-01T00:00:00+08:00"}]},
        id="log",
    ),
    pytest.param(
        upload_api.upload_media,
        ("pid", "eid"),
        {"metrics": [{"key": "img", "index": 0, "data": [], "more": [], "create_time": "2024-01-01T00:00:00+08:00"}]},
        id="media",
    ),
    pytest.param(
        upload_api.upload_columns,
        ("user", "proj"),
        {"columns": {"series": [{"key": "loss"}]}},
        id="column",
    ),
]


@pytest.mark.parametrize("api_func, args, kwargs", _RETRIES_ZERO_CASES)
def test_upload_api_uses_retries_zero(monkeypatch, api_func, args, kwargs):
    """验证所有上传 API 函数都传 retries=0，确保 HTTP 层不自行重试。

    retries=0 是刻意设计：将重试责任上移到 transport 层，
    避免 HTTP 层（urllib3）与 transport 层同时重试造成行为不可预期。
    """
    captured: list = []

    def fake_post(url, data=None, **post_kwargs):
        captured.append(post_kwargs.get("retries"))

    monkeypatch.setattr(upload_api.client, "post", fake_post)

    api_func(*args, **kwargs)

    assert captured == [0]


# ── Section 2: handler → sender.upload 错误传播 ─────────────
# handler 必须让 ApiError 穿透到 sender.upload，后者决定是否重试。
# 任何 handler 不得将 HTTP API 调用包裹在 safe.block 中。


_RETRY_CONTRACT_CASES = [
    pytest.param("scalar", "upload_scalar", _make_scalar_record, id="scalar"),
    pytest.param("log", "upload_log", _make_log_record, id="log"),
    pytest.param("column", "upload_columns", _make_column_record, id="column"),
]


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_5xx_reraises_for_transport_retry(tmp_path, record_type, api_func_name, make_record):
    """5xx ApiError 必须穿透 handler → sender.upload 并 re-raise，供 transport 层重试。

    回归守卫：handler 不得将 HTTP 调用包裹在 safe.block 中，
    否则 sender.upload 永远看不到错误，transport 层无法重试，数据静默丢失。
    """
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(502),
    ):
        with pytest.raises(ApiError):
            sender.upload(record_type, [make_record()])

    assert tracker.snapshot().uploaded_records == 0


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_4xx_swallowed_non_retryable(tmp_path, record_type, api_func_name, make_record):
    """4xx ApiError 被 sender.upload 吞掉（业务拒绝），不重试，records 不递进。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(400),
    ):
        sender.upload(record_type, [make_record()])

    assert tracker.snapshot().uploaded_records == 0


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_success_advances_records(tmp_path, record_type, api_func_name, make_record):
    """成功上传递进 uploaded_records。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)

    with patch(f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}"):
        sender.upload(record_type, [make_record()])

    assert tracker.snapshot().uploaded_records == 1
