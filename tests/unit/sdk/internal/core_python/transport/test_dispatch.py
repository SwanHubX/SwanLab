import threading
from unittest.mock import MagicMock, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.buffer import RecordBuffer
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch


def test_dispatch_groups_by_type(make_scalar_record, make_config_record):
    """混合类型 record 被正确分组分发。"""
    dispatch = Dispatch(cond=threading.Condition(), buf=RecordBuffer())

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(dispatch, "_handle_record_by_type", return_value=(True, [])) as mock_handle:
        dispatch(metric_records + config_records)
        calls = mock_handle.call_args_list
        assert calls[0] == (("metric", metric_records),)
        assert calls[1] == (("config", config_records),)


def test_dispatch_calls_correct_handler(make_scalar_record):
    """各 _handle_{kind} 被调用且参数正确。"""
    dispatch = Dispatch(cond=threading.Condition(), buf=RecordBuffer())
    records = [make_scalar_record(step=1)]

    with patch.object(dispatch, "_handle_record_by_type", return_value=(True, [])) as mock_handle:
        dispatch(records)
        mock_handle.assert_called_once_with("metric", records)


def test_dispatch_error_rollback(make_scalar_record):
    """上传失败回滚到 buffer 头部（原地插入，buffer 对象引用不变）。"""
    buf = RecordBuffer()
    mock_sender = MagicMock()
    mock_sender.upload.side_effect = RuntimeError("upload failed")
    dispatch = Dispatch(cond=threading.Condition(), buf=buf, sender=mock_sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 11

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        dispatch(records)

    drained = buf.drain()
    assert drained == records


def test_dispatch_rollback_skips_records_already_buffered_by_num(make_scalar_record):
    """回滚时若 buffer 中已有同 num record，不再重复插入。"""
    buf = RecordBuffer()
    existing = make_scalar_record(step=1)
    existing.num = 21
    buf.extend([existing])

    dispatch = Dispatch(cond=threading.Condition(), buf=buf)

    duplicate = make_scalar_record(step=1)
    duplicate.num = 21
    other = make_scalar_record(step=2)
    other.num = 22

    with patch.object(dispatch, "_handle_record_by_type", return_value=(False, [duplicate, other])):
        dispatch([duplicate, other])

    assert [record.num for record in buf.drain()] == [22, 21]


def test_dispatch_skips_unknown_type(make_scalar_record):
    """未知 kind 无 handler 时不报错，静默跳过。"""
    dispatch = Dispatch(cond=threading.Condition(), buf=RecordBuffer())

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"unknown_kind": [make_scalar_record(step=1)]})

        with patch("swanlab.sdk.internal.core_python.transport.dispatch.console.warning") as mock_warning:
            dispatch([make_scalar_record(step=1)])
            mock_warning.assert_called_once()


def test_dispatch_mixed_type_partial_failure_rollback(make_scalar_record, make_config_record):
    """混合类型中一种上传失败时，后续类型的 records 也被统一回滚到 buffer。"""
    buf = RecordBuffer()
    mock_sender = MagicMock()

    def upload_side_effect(record_type, records):
        if record_type == "metric":
            raise RuntimeError("upload failed")

    mock_sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(cond=threading.Condition(), buf=buf, sender=mock_sender)

    metric_records = [make_scalar_record(step=1)]
    metric_records[0].num = 31
    config_records = [make_config_record()]
    config_records[0].num = 32

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        dispatch(metric_records + config_records)

    drained = buf.drain()
    assert drained[0] is metric_records[0]
    assert drained[1] is config_records[0]


def test_dispatch_handle_record_type_success_calls_callback(make_scalar_record):
    """上传成功后返回成功状态并触发上传回调。"""
    uploaded = []
    callback = MagicMock()
    sender = MagicMock()

    def upload_side_effect(record_type, records):
        uploaded.append((record_type, list(records)))

    sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(cond=threading.Condition(), buf=RecordBuffer(), upload_callback=callback, sender=sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 41

    success, failed = dispatch._handle_record_by_type("metric", records)

    assert success is True
    assert failed == []
    assert uploaded == [("metric", records)]
    callback.assert_called_once_with(1)


def test_dispatch_handle_record_type_returns_failed_tail_after_retries(make_scalar_record):
    """重试耗尽后返回当前组中尚未成功上传的 records。"""
    sender = MagicMock()
    sender.upload.side_effect = [None, RuntimeError("boom"), RuntimeError("boom"), RuntimeError("boom")]
    dispatch = Dispatch(
        cond=threading.Condition(),
        buf=RecordBuffer(),
        sender=sender,
        max_retries=3,
        initial_backoff=0,
    )

    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 51
    second.num = 52

    with patch(
        "swanlab.sdk.internal.core_python.transport.dispatch.generate_chunks",
        return_value=[([first], 1), ([second], 1)],
    ):
        success, failed = dispatch._handle_record_by_type("metric", [first, second])

    assert success is False
    assert failed == [second]


def test_dispatch_record_types_follow_proto_descriptor():
    """分发类型集合应直接来自 Record.record_type oneof。"""
    expected = frozenset(field.name for field in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)
    assert Dispatch._RECORD_TYPES == expected
