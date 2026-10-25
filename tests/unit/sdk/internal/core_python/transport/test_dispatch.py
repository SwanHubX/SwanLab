from unittest.mock import MagicMock, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord, MediaValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult


def _ok(records) -> UploadResult:
    return UploadResult(confirmed_records=list(records))


def _failed(records, message: str = "upload failed") -> UploadResult:
    return UploadResult(
        remaining_records=list(records),
        failures=[UploadFailure(message=message, records=list(records))],
    )


def _make_media_record(filename: str = "image.png") -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        media=MediaRecord(
            key="examples/image",
            step=1,
            type=ColumnType.COLUMN_TYPE_IMAGE,
            timestamp=timestamp,
            value=MediaValue(items=[MediaItem(filename=filename)]),
        )
    )


def test_dispatch_groups_by_type(mock_ctx, make_scalar_record, make_config_record):
    """混合类型 record 按类型分组分发。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(
        dispatch, "_upload_record_type", side_effect=[_ok(metric_records), _ok(config_records)]
    ) as mock_handle:
        result = dispatch(metric_records + config_records)
        calls = mock_handle.call_args_list
        assert result.remaining_records == []
        assert result.failures == []
        assert result.confirmed_records == metric_records + config_records
        assert calls[0] == (("scalar", metric_records),)
        assert calls[1] == (("save", config_records),)


@pytest.mark.parametrize(
    "error, terminal",
    [
        (RuntimeError("upload exploded"), False),
        (OSError("disk unavailable"), True),
    ],
    ids=["runtime-error", "os-error"],
)
def test_dispatch_sender_unexpected_error_converted_to_failure(mock_ctx, make_scalar_record, error, terminal):
    """sender 未预期异常转为失败结果，本地错误终结，内部错误保留重试。"""
    mock_sender = MagicMock()
    mock_sender.upload.side_effect = error
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=mock_sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 11

    result = dispatch(records)

    if terminal:
        assert result.confirmed_records == records
        assert result.remaining_records == []
    else:
        assert result.confirmed_records == []
        assert result.remaining_records == records
    failure = result.failures[0]
    assert type(error).__name__ in failure.message
    assert failure.exception is error


def test_dispatch_confirms_untypeable_records(mock_ctx):
    """不可识别记录跳过上传，但已终止处理，计入 confirmed 与 uploaded_records。"""
    tracker = UploadTracker()
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, tracker=tracker)
    records = [Record(), Record()]

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        result = dispatch(records)

    assert result.confirmed_records == records
    assert result.remaining_records == []
    assert result.failures == []
    assert tracker.snapshot().uploaded_records == 2


def test_dispatch_settles_non_retryable_sender_error(mock_ctx, make_scalar_record):
    """sender 抛确定性协议错误时终结当前分片，不交给 worker 退避重试。"""
    tracker = UploadTracker()
    sender = MagicMock()
    sender.upload.side_effect = ApiError(
        code="invalid_presigned_response",
        message="presigned response returned 1 urls for 2 paths",
        retryable=False,
    )
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender, tracker=tracker)

    result = dispatch(records)

    assert result.confirmed_records == records
    assert result.remaining_records == []
    # confirmed 表示“已终止处理”，确定性失败同样计入
    assert tracker.snapshot().uploaded_records == 2


def test_dispatch_advances_uploaded_for_confirmed_records(mock_ctx, make_scalar_record):
    """uploaded_records 统计已终止处理数，含终结失败的记录。"""
    tracker = UploadTracker()
    records = [make_scalar_record(step=i) for i in range(3)]
    for index, record in enumerate(records):
        record.num = 300 + index

    def _upload(record_type, chunk):
        return UploadResult(
            confirmed_records=list(chunk),
            failures=[UploadFailure(message="http 400", records=[chunk[0]])],
        )

    sender = MagicMock()
    sender.upload.side_effect = _upload
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender, tracker=tracker)

    result = dispatch(records)

    assert len(result.confirmed_records) == 3
    assert tracker.snapshot().uploaded_records == 3


def test_dispatch_settled_type_does_not_block_next_type(mock_ctx, make_scalar_record):
    """某个类型确定性失败后，同批后续类型仍继续上传，不被失败文件阻塞。"""
    media = _make_media_record()
    scalar = make_scalar_record(step=1)

    def _upload(record_type, chunk):
        if record_type == "media":
            return UploadResult(
                confirmed_records=list(chunk),
                failures=[UploadFailure(message="presigned url mismatch", records=list(chunk))],
            )
        return _ok(chunk)

    sender = MagicMock()
    sender.upload.side_effect = _upload
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender)

    result = dispatch([media, scalar])

    assert result.confirmed_records == [media, scalar]
    assert result.remaining_records == []


@pytest.mark.parametrize("position", [0, 1, 2], ids=["head", "middle", "tail"])
def test_dispatch_keeps_valid_records_when_untypeable(mock_ctx, make_scalar_record, position):
    """一条不可识别 record 只跳过自身，有效记录继续上传。"""
    valid = [make_scalar_record(step=1), make_scalar_record(step=2)]
    valid[0].num, valid[1].num = 101, 102
    batch: list[Record] = list(valid)
    batch.insert(position, Record())

    sender = MagicMock()
    sender.upload.return_value = _ok(valid)
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender)

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        result = dispatch(batch)

    assert sender.upload.call_count == 1
    assert sender.upload.call_args.args[1] == valid
    assert result.confirmed_records == [batch[position]] + valid
    assert result.remaining_records == []
    assert result.failures == []


def test_dispatch_propagates_type_error(mock_ctx, make_scalar_record):
    """非 Record 属于调用方 bug，直接抛错，不做整批恢复。"""
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch)
    record = make_scalar_record(step=1)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type",
            side_effect=TypeError("Expected Record instance"),
        ),
        pytest.raises(TypeError, match="Expected Record instance"),
    ):
        dispatch([record])


def test_dispatch_mixed_type_partial_failure_rollback(mock_ctx, make_scalar_record, make_config_record):
    """混合类型中前序组失败时保留当前组与后续组。"""
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

    assert result.remaining_records
    assert result.remaining_records[0] is metric_records[0]
    assert result.remaining_records[1] is config_records[0]
    assert result.failures[0].message == "upload failed"


def test_dispatch_advances_tracker_by_confirmed_records(mock_ctx, make_scalar_record):
    """根据 confirmed_records 推进 tracker 统计。"""
    tracker = UploadTracker()
    sender = MagicMock()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 41
    records[1].num = 42
    sender.upload.return_value = _ok(records)
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender, tracker=tracker)

    result = dispatch(records)

    assert result.remaining_records == []
    assert result.failures == []
    assert tracker.snapshot().uploaded_records == 2


def test_dispatch_partial_chunk_success_advances_tracker_and_keeps_tail(mock_ctx, make_scalar_record):
    """前序分片成功、后续分片失败时，确认成功分片并保留剩余分片。"""
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


def test_dispatch_first_chunk_success_second_chunk_exception_keeps_prefix(mock_ctx, make_scalar_record):
    """前序分片成功、后续分片抛异常时，确认成功前缀并保留剩余分片。"""
    tracker = UploadTracker()
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    third = make_scalar_record(step=3)
    first.num, second.num, third.num = 51, 52, 53
    sender.upload.side_effect = [_ok([first]), RuntimeError("chunk exploded")]
    dispatch = Dispatch(batch_size=1, sender=sender, tracker=tracker)

    result = dispatch([first, second, third])

    assert result.confirmed_records == [first]
    assert result.remaining_records == [second, third]
    failure = result.failures[0]
    assert failure.exception is not None
    assert tracker.snapshot().uploaded_records == 1


def test_dispatch_scalar_success_media_exception_keeps_confirmed_prefix(mock_ctx, make_scalar_record):
    """scalar 分组成功、media 分组抛异常时，确认 scalar 并保留 media。"""
    tracker = UploadTracker()
    sender = MagicMock()

    def upload(record_type, records):
        if record_type == "scalar":
            return _ok(records)
        raise RuntimeError("media upload exploded")

    sender.upload.side_effect = upload
    dispatch = Dispatch(batch_size=mock_ctx.config.record_batch, sender=sender, tracker=tracker)
    scalar = make_scalar_record(step=1)
    scalar.num = 61
    media = _make_media_record()
    media.num = 62

    result = dispatch([scalar, media])

    assert result.confirmed_records == [scalar]
    assert result.remaining_records == [media]
    assert tracker.snapshot().uploaded_records == 1


def test_dispatch_next_round_retries_only_remaining(mock_ctx, make_scalar_record):
    """下一轮重试剩余记录，不重放成功前缀。"""
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num, second.num = 71, 72
    sender.upload.side_effect = [_ok([first]), RuntimeError("chunk exploded"), _ok([second])]
    dispatch = Dispatch(batch_size=1, sender=sender)

    first_round = dispatch([first, second])
    assert first_round.confirmed_records == [first]
    assert first_round.remaining_records == [second]

    second_round = dispatch(first_round.remaining_records)

    assert second_round.remaining_records == []
    assert second_round.failures == []
    uploaded_chunks = [call.args[1] for call in sender.upload.call_args_list]
    assert uploaded_chunks == [[first], [second], [second]]


def test_dispatch_result_contract_violation_retains_chunk(mock_ctx, make_scalar_record):
    """sender 返回结果记录数不匹配时保留当前及后续分片。"""
    tracker = UploadTracker()
    sender = MagicMock()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num, second.num = 61, 62
    sender.upload.return_value = UploadResult()
    dispatch = Dispatch(batch_size=1, sender=sender, tracker=tracker)

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        result = dispatch([first, second])

    assert result.remaining_records == [first, second]
    assert tracker.snapshot().uploaded_records == 0
