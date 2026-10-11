"""
重试设计契约测试：验证上传链路中 ApiError 的处理规则和 retries=0 约定。

上传错误处理分四层：
  1. API 层（api/upload.py）：所有 /house/metrics 和 column 上传使用 retries=0，
     HTTP 层不重试，重试交由 transport 层处理。
  2. sender.upload：错误决定点。ApiError / requests 异常转为结果内的 UploadFailure，
     记录保留在 remaining，不向上抛出、不吞错、不消费。
  3. Dispatch：汇总 UploadResult，依据 confirmed_records 推进统计，不吞错误。
  4. worker：持有未确认记录，节流告警，按重试策略再次调度（P3）。

契约：handler 不得将 HTTP API 调用包裹在 safe.block 中，
否则 sender.upload 拿不到异常详情，失败原因无法结构化。
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
from swanlab.sdk.internal.core_python.api import save as save_api
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
    """抑制 sender 与 safe.block 的 console 输出。"""
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
    pytest.param(
        upload_api.upload_conda,
        ("user", "proj", "eid"),
        {"content": "conda-list"},
        id="conda",
    ),
    pytest.param(
        upload_api.upload_requirements,
        ("user", "proj", "eid"),
        {"content": "numpy==1.0"},
        id="requirements",
    ),
    pytest.param(
        upload_api.upload_metadata,
        ("user", "proj", "eid"),
        {"content": {"hostname": "gpu-01"}},
        id="metadata",
    ),
    pytest.param(
        upload_api.upload_config,
        ("user", "proj", "eid"),
        {"content": {"lr": {"value": 0.01}}},
        id="config",
    ),
    pytest.param(
        save_api.prepare_save_files,
        ("eid",),
        {"files": [{"path": "model.bin", "size": 8, "md5": "0", "mimeType": "text/plain"}]},
        id="save-prepare",
    ),
    pytest.param(
        save_api.complete_save_files,
        ("eid",),
        {"files": [{"path": "model.bin", "state": "UPLOADED"}]},
        id="save-complete",
    ),
    pytest.param(
        save_api.prepare_multipart_save,
        ("eid",),
        {"files": [{"path": "large.bin", "size": 8, "md5": "0", "mimeType": "text/plain", "count": 1}]},
        id="save-prepare-multipart",
    ),
    pytest.param(
        save_api.complete_multipart_save,
        ("eid",),
        {"files": [{"path": "large.bin", "parts": [], "uploadId": "uid"}]},
        id="save-complete-multipart",
    ),
]


@pytest.mark.parametrize("api_func, args, kwargs", _RETRIES_ZERO_CASES)
def test_upload_api_uses_retries_zero(monkeypatch, api_func, args, kwargs):
    """验证所有上传 API 函数都传 retries=0，确保 HTTP 层不自行重试。

    retries=0 是刻意设计：将重试责任上移到 transport 层的退避调度，
    避免 HTTP 层（urllib3）与 transport 层同时重试造成行为不可预期。
    """
    captured: list = []

    class _FakeResponse:
        data = {}

    def fake_request(url, data=None, **request_kwargs):
        captured.append(request_kwargs.get("retries"))
        return _FakeResponse()

    monkeypatch.setattr(upload_api.client, "post", fake_request)
    monkeypatch.setattr(upload_api.client, "put", fake_request)

    api_func(*args, **kwargs)

    assert captured == [0]


def test_presigned_put_uses_retries_zero(monkeypatch):
    """upload_resource 的预签名地址请求同样单次执行。"""

    class _FakePresignResponse:
        data = {"urls": []}

    captured: list = []

    def fake_post(url, data=None, **post_kwargs):
        captured.append(post_kwargs.get("retries"))
        return _FakePresignResponse()

    monkeypatch.setattr(upload_api.client, "post", fake_post)
    upload_api.upload_resource(MagicMock(), "eid", paths=[], buffers=[])

    assert captured == [0]


# ── Section 2: handler → sender.upload 错误转结果 ─────────────
# handler 必须让异常穿透到 sender.upload，由后者转为结果内的 UploadFailure。
# 任何 handler 不得将 HTTP API 调用包裹在 safe.block 中，
# 否则 sender.upload 拿不到异常详情，失败原因无法结构化、记录无法如实保留。


_RETRY_CONTRACT_CASES = [
    pytest.param("scalar", "upload_scalar", _make_scalar_record, id="scalar"),
    pytest.param("log", "upload_log", _make_log_record, id="log"),
    pytest.param("column", "upload_columns", _make_column_record, id="column"),
]


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_5xx_becomes_failure_result(tmp_path, record_type, api_func_name, make_record):
    """5xx ApiError 转为结果内的失败：记录全部保留在 remaining，不向上抛出。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = make_record()

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(502),
    ):
        result = sender.upload(record_type, [record])

    assert result.confirmed_records == []
    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.status == 502
    assert failure.source == "swanlab_api"
    assert failure.stage == "submit_metrics"
    assert failure.exception is not None
    # sender 不推进统计，统计由 Dispatch 依据 confirmed_records 推进
    assert tracker.snapshot().uploaded_records == 0


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_4xx_becomes_failure_result(tmp_path, record_type, api_func_name, make_record):
    """4xx 业务拒绝同样保留记录并报告失败：不吞错、不消费、不冒充成功。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = make_record()

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(400),
    ):
        result = sender.upload(record_type, [record])

    assert result.confirmed_records == []
    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.status == 400
    assert failure.source == "swanlab_api"
    assert tracker.snapshot().uploaded_records == 0


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_success_returns_confirmed(tmp_path, record_type, api_func_name, make_record):
    """成功上传返回 confirmed，不抛异常；统计推进由 Dispatch 负责。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = make_record()

    with patch(f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}"):
        result = sender.upload(record_type, [record])

    assert result.ok
    assert result.confirmed_records == [record]
    assert result.remaining_records == []
    assert tracker.snapshot().uploaded_records == 0
