from typing import Sequence, cast
from unittest.mock import MagicMock, call, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.sender import (
    HttpRecordTransport,
    create_record_transport,
    trace_records,
    upload_records,
)


def test_http_record_transport_upload_record_group_skips_empty():
    """验证上传空列表时不抛异常，直接跳过。"""
    transport = HttpRecordTransport()
    transport.upload_record_group("metric", [])  # should not raise


def test_http_record_transport_close():
    """验证 close() 不抛异常（当前为空实现）。"""
    transport = HttpRecordTransport()
    transport.close()  # should not raise


# ─────────────────── create_record_transport ───────────────────


def test_create_record_transport_returns_http_transport():
    """验证工厂函数返回 HttpRecordTransport 实例。"""
    transport = create_record_transport()

    assert isinstance(transport, HttpRecordTransport)


# ─────────────────── trace_records ───────────────────


def test_trace_records_groups_records_by_record_type(make_scalar_record, make_config_record):
    """验证 trace_records 按 record_type 分组后分别上传，同类型 Record 合并为一次调用。"""
    transport = MagicMock()
    metric_1 = make_scalar_record(step=1)
    config = make_config_record()
    metric_2 = make_scalar_record(step=2)
    records = [metric_1, config, metric_2]

    trace_records(records, per_request_len=-1, transport=transport)

    assert transport.upload_record_group.call_args_list == [
        call("metric", [metric_1, metric_2]),
        call("config", [config]),
    ]
    transport.close.assert_not_called()


def test_trace_records_uploads_all_records_in_batches(make_scalar_record):
    """验证超过 per_request_len 的记录被正确分片，每片独立上传并触发回调。"""
    transport = MagicMock()
    callback = MagicMock()
    records = [make_scalar_record(step=index) for index in range(2500)]

    with patch(
        "swanlab.sdk.internal.core_python.uploader.sender.create_record_transport", return_value=transport
    ) as factory:
        trace_records(records, per_request_len=1000, upload_callback=callback)

    assert transport.upload_record_group.call_args_list == [
        call("metric", records[:1000]),
        call("metric", records[1000:2000]),
        call("metric", records[2000:]),
    ]
    assert callback.call_args_list == [call(1000), call(1000), call(500)]
    factory.assert_called_once_with()
    transport.close.assert_called_once_with()


def test_trace_records_rejects_serialized_bytes(make_scalar_record):
    """验证传入序列化字节而非 Record 实例时抛出 TypeError。"""
    transport = MagicMock()
    bad_records = cast(Sequence[Record], [make_scalar_record(step=7).SerializeToString()])

    try:
        trace_records(bad_records, transport=transport)
    except TypeError as exc:
        assert "Record" in str(exc)
    else:
        raise AssertionError("trace_records should reject serialized bytes")


def test_trace_records_none_input():
    """验证 records=None 时直接返回，不触发任何上传。"""
    transport = MagicMock()
    trace_records(records=None, transport=transport)

    transport.upload_record_group.assert_not_called()


def test_trace_records_empty_input():
    """验证空列表时直接返回，不触发任何上传。"""
    transport = MagicMock()
    trace_records(records=[], transport=transport)

    transport.upload_record_group.assert_not_called()


def test_trace_records_creates_and_closes_transport_when_none_provided(make_scalar_record):
    """验证未提供 transport 时自动创建并在 finally 中关闭。"""
    records = [make_scalar_record(step=1)]

    with patch(
        "swanlab.sdk.internal.core_python.uploader.sender.create_record_transport"
    ) as factory:
        mock_transport = MagicMock()
        factory.return_value = mock_transport

        trace_records(records, per_request_len=-1)

    factory.assert_called_once_with()
    mock_transport.upload_record_group.assert_called_once()
    mock_transport.close.assert_called_once_with()


def test_trace_records_does_not_close_provided_transport(make_scalar_record):
    """验证外部提供的 transport 在完成后不会被关闭。"""
    transport = MagicMock()
    records = [make_scalar_record(step=1)]

    trace_records(records, per_request_len=-1, transport=transport)

    transport.upload_record_group.assert_called_once()
    transport.close.assert_not_called()


def test_trace_records_callback_none(make_scalar_record):
    """验证 upload_callback=None 时不触发回调。"""
    transport = MagicMock()
    records = [make_scalar_record(step=1)]

    trace_records(records, per_request_len=-1, upload_callback=None, transport=transport)

    transport.upload_record_group.assert_called_once()


# ─────────────────── upload_records ───────────────────


def test_upload_records_empty_list():
    """验证 upload_records 传入空列表时输出 debug 日志并返回。"""
    with patch("swanlab.sdk.internal.core_python.uploader.sender.console") as mock_console:
        upload_records([])

    mock_console.debug.assert_called_once_with("No records to upload.")


def test_upload_records_delegates_to_trace_records(make_scalar_record):
    """验证 upload_records 正确委托给 trace_records。"""
    records = [make_scalar_record(step=1)]
    callback = MagicMock()

    with patch("swanlab.sdk.internal.core_python.uploader.sender.trace_records") as mock_trace:
        upload_records(records, upload_callback=callback, per_request_len=500)

    mock_trace.assert_called_once_with(records, per_request_len=500, upload_callback=callback)
