"""上传 API retries=0 约定与 sender 错误转换测试。"""

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


class _FakeApiErrorResponse:
    """ApiError 测试替身。"""

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
    """抑制测试期间的 console 输出。"""
    with (
        patch("swanlab.sdk.internal.pkg.console.trace"),
        patch("swanlab.sdk.internal.pkg.console.warning"),
    ):
        yield


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
    """上传 API 函数传入 retries=0，由 transport 层控制重试。"""
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




_RETRY_CONTRACT_CASES = [
    pytest.param("scalar", "upload_scalar", _make_scalar_record, id="scalar"),
    pytest.param("log", "upload_log", _make_log_record, id="log"),
    pytest.param("column", "upload_columns", _make_column_record, id="column"),
]


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
@pytest.mark.parametrize("status", [400, 502], ids=["4xx", "5xx"])
def test_http_error_becomes_failure_result(tmp_path, record_type, api_func_name, make_record, status):
    """HTTP 失败时保留记录与错误详情，不推进 tracker 统计。"""
    tracker = UploadTracker()
    sender = _make_sender(tmp_path)
    sender.set_tracker(tracker)
    record = make_record()

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(status),
    ):
        result = sender.upload(record_type, [record])

    assert result.confirmed_records == []
    assert result.remaining_records == [record]
    failure = result.failures[0]
    assert failure.status == status
    assert failure.source == "swanlab_api"
    assert failure.stage == "submit_metrics"
    assert failure.code == "test-error"
    assert failure.exception is not None
    assert tracker.snapshot().uploaded_records == 0


@pytest.mark.parametrize("record_type, api_func_name, make_record", _RETRY_CONTRACT_CASES)
def test_success_returns_confirmed(tmp_path, record_type, api_func_name, make_record):
    """上传成功时返回 confirmed_records，统计推进由 Dispatch 负责。"""
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


@pytest.mark.parametrize(
    "record_type, api_func_name, make_good_record, make_bad_record",
    [
        pytest.param("scalar", "upload_scalar", _make_scalar_record, lambda: Record(), id="scalar"),
        pytest.param("log", "upload_log", _make_log_record, lambda: Record(), id="log"),
    ],
)
def test_mixed_batch_http_failure_splits_failure_sources(
    tmp_path, record_type, api_func_name, make_good_record, make_bad_record
):
    """混合批次遇到 HTTP 失败时，有效记录保留重试，无效记录视为已终结不再重试，失败来源拆分。"""
    sender = _make_sender(tmp_path)
    good = make_good_record()
    bad = make_bad_record()

    with patch(
        f"swanlab.sdk.internal.core_python.transport.sender.{api_func_name}",
        side_effect=_make_api_error(503),
    ):
        result = sender.upload(record_type, [good, bad])

    assert result.confirmed_records == []
    # 只有有效记录保留重试；无效记录（blocked）不进 remaining
    assert result.remaining_records == [good]
    assert bad not in result.remaining_records

    sources = {f.source for f in result.failures}
    assert "swanlab_api" in sources
    assert "local_data" in sources

    api_failure = next(f for f in result.failures if f.source == "swanlab_api")
    local_failure = next(f for f in result.failures if f.source == "local_data")

    assert good in api_failure.records
    assert bad not in api_failure.records
    assert bad in local_failure.records
    assert good not in local_failure.records
