from unittest.mock import MagicMock, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult


def _ok(records) -> UploadResult:
    return UploadResult(confirmed_records=list(records))


def _failed(records, message: str = "upload failed") -> UploadResult:
    return UploadResult(
        remaining_records=list(records),
        failures=[UploadFailure(source="swanlab_api", message=message, records=list(records))],
    )


def test_dispatch_groups_by_type(mock_ctx, make_scalar_record, make_config_record):
    """混合类型 record 被正确分组分发。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(
        dispatch, "_upload_record_type", side_effect=[_ok(metric_records), _ok(config_records)]
    ) as mock_handle:
        result = dispatch(metric_records + config_records)
        calls = mock_handle.call_args_list
        assert result.ok
        assert result.confirmed_records == metric_records + config_records
        assert calls[0] == (("scalar", metric_records),)
        assert calls[1] == (("save", config_records),)


def test_dispatch_calls_correct_handler(mock_ctx, make_scalar_record):
    """各 record_type 被正确交给 _upload_record_type()。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)
    records = [make_scalar_record(step=1)]

    with patch.object(dispatch, "_upload_record_type", return_value=_ok(records)) as mock_handle:
        dispatch(records)
        mock_handle.assert_called_once_with("scalar", records)


def test_dispatch_sender_unexpected_error_propagates(mock_ctx, make_scalar_record):
    """sender 的未预期异常（编程错误）不被 Dispatch 吞掉，向上交给 worker 兜底。"""
    mock_sender = MagicMock()
    mock_sender.upload.side_effect = RuntimeError("upload failed")
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=mock_sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 11

    try:
        dispatch(records)
    except RuntimeError:
        pass
    else:
        raise AssertionError("RuntimeError should propagate out of Dispatch")


def test_dispatch_failure_tail_keeps_failed_and_unprocessed_order(mock_ctx, make_scalar_record, make_config_record):
    """当前组失败部分和后续未处理组保持原顺序返回。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)

    failed_metric = make_scalar_record(step=1)
    failed_metric.num = 21
    later_config = make_config_record()
    later_config.num = 22

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"metric": [failed_metric], "config": [later_config]})

        with patch.object(dispatch, "_upload_record_type", return_value=_failed([failed_metric])):
            result = dispatch([failed_metric, later_config])

    assert not result.ok
    assert result.remaining_records == [failed_metric, later_config]
    assert result.confirmed_records == []


def test_dispatch_retains_unknown_type(mock_ctx, make_scalar_record):
    """未知 kind 无 handler：整组保留并报告失败，不静默跳过冒充完成。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)
    record = make_scalar_record(step=1)

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"unknown_kind": [record]})

        with patch("swanlab.sdk.internal.core_python.transport.dispatch.console.warning") as mock_warning:
            result = dispatch([record])

    assert result.remaining_records == [record]
    assert result.failures[0].source == "local_data"
    assert "unknown record type" in result.failures[0].message
    mock_warning.assert_called_once()


def test_dispatch_retains_untypeable_records(mock_ctx):
    """无法确定 record_type 的记录：整批保留并报告，不静默丢弃。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)
    record = Record()  # 未设置任何 record_type

    result = dispatch([record])

    assert result.remaining_records == [record]
    assert len(result.failures) == 1
    assert result.failures[0].source == "local_data"


def test_dispatch_mixed_type_partial_failure_rollback(mock_ctx, make_scalar_record, make_config_record):
    """混合类型中一种上传失败时，返回当前失败组和后续组，成功部分不回滚。"""
    mock_sender = MagicMock()

    def upload_side_effect(record_type, records):
        if record_type == "scalar":
            return _failed(records)
        return _ok(records)

    mock_sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=mock_sender)

    metric_records = [make_scalar_record(step=1)]
    metric_records[0].num = 31
    config_records = [make_config_record()]
    config_records[0].num = 32

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        result = dispatch(metric_records + config_records)

    assert not result.ok
    assert result.remaining_records[0] is metric_records[0]
    assert result.remaining_records[1] is config_records[0]
    assert result.failures[0].message == "upload failed"


def test_dispatch_advances_tracker_by_confirmed_records(mock_ctx, make_scalar_record):
    """统计由 Dispatch 依据 confirmed_records 推进：成功确认的数量恰好推进一次。"""
    tracker = UploadTracker()
    sender = MagicMock()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 41
    records[1].num = 42
    sender.upload.return_value = _ok(records)
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender, tracker=tracker)

    result = dispatch(records)

    assert result.ok
    assert tracker.snapshot().uploaded_records == 2


def test_dispatch_partial_chunk_success_advances_tracker_and_keeps_tail(mock_ctx, make_scalar_record):
    """前一个 chunk 成功、后一个失败：成功部分推进统计，失败 chunk 及后续保留。"""
    tracker = UploadTracker()
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    third = make_scalar_record(step=3)
    first.num, second.num, third.num = 51, 52, 53
    sender.upload.side_effect = [_ok([first]), _failed([second])]
    dispatch = Dispatch(batch_size=1, sender=sender, tracker=tracker)

    result = dispatch([first, second, third])

    assert result.confirmed_records == [first]
    assert result.remaining_records == [second, third]
    assert tracker.snapshot().uploaded_records == 1


def test_dispatch_handle_record_type_returns_failed_tail_on_chunk_failure(mock_ctx, make_scalar_record):
    """某个 chunk 上传失败时返回当前组中尚未成功上传的 records。"""
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 51
    second.num = 52
    sender.upload.side_effect = [_ok([first]), _failed([second])]
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender)

    with patch(
        "swanlab.sdk.internal.core_python.transport.dispatch.generate_chunks",
        return_value=[([first], 1), ([second], 1)],
    ):
        result = dispatch._upload_record_type("scalar", [first, second])

    assert result.confirmed_records == [first]
    assert result.remaining_records == [second]


def test_dispatch_result_contract_violation_retains_chunk(mock_ctx, make_scalar_record):
    """sender 结果未覆盖全部输入（契约违约）：整段保留并报告，不静默丢数据。"""
    tracker = UploadTracker()
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num, second.num = 61, 62
    # 返回空结果：1 条输入没有任何归属，违反结果契约
    sender.upload.return_value = UploadResult()
    dispatch = Dispatch(batch_size=1, sender=sender, tracker=tracker)

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        result = dispatch([first, second])

    assert result.remaining_records == [first, second]
    assert result.failures[0].source == "worker_internal"
    assert tracker.snapshot().uploaded_records == 0


def test_dispatch_record_types_follow_proto_descriptor():
    """分发类型集合应直接来自 Record.record_type oneof。"""
    expected = frozenset(field.name for field in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)
    assert Dispatch._RECORD_TYPES == expected
