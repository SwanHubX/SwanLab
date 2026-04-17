import threading
from unittest.mock import MagicMock, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch


def test_dispatch_groups_by_type(make_scalar_record, make_config_record):
    """混合类型 record 被正确分组分发。"""
    dispatch = Dispatch()

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(dispatch, "_upload_record_type", return_value=(True, [])) as mock_handle:
        success, failed = dispatch(metric_records + config_records)
        calls = mock_handle.call_args_list
        assert success is True
        assert failed == []
        assert calls[0] == (("metric", metric_records),)
        assert calls[1] == (("config", config_records),)


def test_dispatch_calls_correct_handler(make_scalar_record):
    """各 record_type 被正确交给 _upload_record_type()。"""
    dispatch = Dispatch()
    records = [make_scalar_record(step=1)]

    with patch.object(dispatch, "_upload_record_type", return_value=(True, [])) as mock_handle:
        dispatch(records)
        mock_handle.assert_called_once_with("metric", records)


def test_dispatch_returns_failed_records_without_mutating_external_buffer(make_scalar_record):
    """上传失败时返回待重试 records，由上层决定如何保留。"""
    mock_sender = MagicMock()
    mock_sender.upload.side_effect = RuntimeError("upload failed")
    dispatch = Dispatch(sender=mock_sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 11

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        success, failed = dispatch(records)

    assert success is False
    assert failed == records


def test_dispatch_failure_tail_keeps_failed_and_unprocessed_order(make_scalar_record, make_config_record):
    """当前组失败部分和后续未处理组保持原顺序返回。"""
    dispatch = Dispatch()

    failed_metric = make_scalar_record(step=1)
    failed_metric.num = 21
    later_config = make_config_record()
    later_config.num = 22

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"metric": [failed_metric], "config": [later_config]})

        with patch.object(dispatch, "_upload_record_type", side_effect=[(False, [failed_metric])]):
            success, failed = dispatch([failed_metric, later_config])

    assert success is False
    assert failed == [failed_metric, later_config]


def test_dispatch_skips_unknown_type(make_scalar_record):
    """未知 kind 无 handler 时不报错，静默跳过。"""
    dispatch = Dispatch()

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"unknown_kind": [make_scalar_record(step=1)]})

        with patch("swanlab.sdk.internal.core_python.transport.dispatch.console.warning") as mock_warning:
            dispatch([make_scalar_record(step=1)])
            mock_warning.assert_called_once()


def test_dispatch_mixed_type_partial_failure_rollback(make_scalar_record, make_config_record):
    """混合类型中一种上传失败时，返回当前失败组和后续组。"""
    mock_sender = MagicMock()

    def upload_side_effect(record_type, records):
        if record_type == "metric":
            raise RuntimeError("upload failed")

    mock_sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(sender=mock_sender)

    metric_records = [make_scalar_record(step=1)]
    metric_records[0].num = 31
    config_records = [make_config_record()]
    config_records[0].num = 32

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        success, failed = dispatch(metric_records + config_records)

    assert success is False
    assert failed[0] is metric_records[0]
    assert failed[1] is config_records[0]


def test_dispatch_handle_record_type_success_calls_callback(make_scalar_record):
    """上传成功后返回成功状态并触发上传回调。"""
    uploaded = []
    callback = MagicMock()
    sender = MagicMock()

    def upload_side_effect(record_type, records):
        uploaded.append((record_type, list(records)))

    sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(upload_callback=callback, sender=sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 41

    success, failed = dispatch._upload_record_type("metric", records)

    assert success is True
    assert failed == []
    assert uploaded == [("metric", records)]
    callback.assert_called_once_with(1)


def test_dispatch_handle_record_type_returns_failed_tail_after_retries(make_scalar_record):
    """重试耗尽后返回当前组中尚未成功上传的 records。"""
    sender = MagicMock()
    sender.upload.side_effect = [None, RuntimeError("boom"), RuntimeError("boom"), RuntimeError("boom")]
    dispatch = Dispatch(
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
        success, failed = dispatch._upload_record_type("metric", [first, second])

    assert success is False
    assert failed == [second]


def test_dispatch_record_types_follow_proto_descriptor():
    """分发类型集合应直接来自 Record.record_type oneof。"""
    expected = frozenset(field.name for field in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)
    assert Dispatch._RECORD_TYPES == expected
