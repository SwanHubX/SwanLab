"""BackgroundConsumer 上行失败分类：可重试回塞、终态丢弃、有界退避。

关键契约：
- 可重试错误必须回塞原批次，退避到期后重新提交，记录不丢；
- 终态错误只告警一次并终结该批，不回塞；
- 单次 flush 按 batch_size 分片提交，失败只回塞未成功的部分；
- scalar / media / histogram 各自独立通道，任一通道失败只影响自己那批；
- 退避通过 queue.get(timeout=...) 表达，不 sleep 轮询。
"""

import io
import os
import queue
import sys
import threading
from typing import cast
from unittest.mock import MagicMock

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord, MediaRecord, ScalarRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel
from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus.events import LogEvent, MetricLogEvent
from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.pkg import console as pkg_console
from swanlab.sdk.internal.run.components.terminal.capture import StreamCapture
from swanlab.sdk.protocol import CoreProtocolError


def _scalar(key: str = "loss") -> ScalarRecord:
    return ScalarRecord(key=key, step=1)


def _histogram(key: str = "weights") -> HistogramRecord:
    return HistogramRecord(key=key, step=1)


def test_retryable_error_restores_batch_and_schedules_backoff(make_retry_consumer):
    consumer, core, clock = make_retry_consumer()
    record = _scalar()
    consumer._scalar_batch.append(record)
    core.upsert_scalars.side_effect = [CoreProtocolError("busy", retryable=True), None]

    consumer._flush()

    assert consumer._scalar_batch == [record]
    assert consumer._retry_at == clock.now + 0.5
    assert consumer._retry_due() is False

    clock.advance(0.5)
    assert consumer._retry_due() is True
    consumer._flush()

    assert consumer._scalar_batch == []
    assert core.upsert_scalars.call_count == 2
    assert core.upsert_scalars.call_args_list[0].args[0] == [record]
    assert core.upsert_scalars.call_args_list[1].args[0] == [record]
    assert consumer._retry_at is None
    assert consumer._retry_attempt == 0


def test_retryable_error_logs_no_warning(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer()
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("busy", retryable=True)

    consumer._flush()

    assert warnings == []


def test_terminal_error_drops_batch_with_one_warning(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, clock = make_retry_consumer()
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("bad record", retryable=False)

    consumer._flush()

    assert consumer._scalar_batch == []
    assert consumer._retry_at is None
    assert len(warnings) == 1

    consumer._flush()  # 空批次不再产生日志
    assert len(warnings) == 1


def test_unknown_exception_is_terminal_and_dropped(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer()
    record = _scalar()
    consumer._scalar_batch.append(record)
    core.upsert_scalars.side_effect = RuntimeError("socket closed")

    consumer._flush()

    assert consumer._scalar_batch == []
    assert consumer._retry_at is None
    assert len(warnings) == 1


def test_buffer_cap_drops_new_records_and_rate_limits_warning(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, _, clock = make_retry_consumer(max_buffered_records=2)

    for index in range(5):
        consumer._append_data_record(_scalar(key=f"loss/{index}"))

    assert [record.key for record in consumer._scalar_batch] == ["loss/0", "loss/1"]
    assert len(warnings) == 1
    assert "Dropping 1 buffered record(s)" in warnings[0][0]

    clock.advance(10.0)
    consumer._append_data_record(_scalar(key="loss/overflow"))
    assert len(warnings) == 2


def test_buffer_byte_cap_rejects_single_oversized_record(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, _, _ = make_retry_consumer(max_buffered_bytes=1)

    consumer._append_data_record(_scalar())

    assert consumer._scalar_batch == []
    assert consumer._buffered_bytes == 0
    assert len(warnings) == 1


def test_buffer_cap_resets_after_flush(make_retry_consumer):
    consumer, core, _ = make_retry_consumer(max_buffered_records=2)
    for index in range(3):
        consumer._append_data_record(_scalar(key=f"loss/{index}"))
    assert len(consumer._scalar_batch) == 2

    consumer._flush()
    assert core.upsert_scalars.call_count == 1
    assert consumer._buffered_bytes == 0

    consumer._append_data_record(_scalar(key="loss/after"))
    assert [record.key for record in consumer._scalar_batch] == ["loss/after"]


def test_terminal_group_does_not_block_other_groups(monkeypatch, make_retry_consumer):
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: None)
    consumer, core, _ = make_retry_consumer()
    consumer._scalar_batch.append(_scalar())
    consumer._log_batch.append(object())  # type: ignore[arg-type]
    core.upsert_scalars.side_effect = CoreProtocolError("bad record", retryable=False)

    consumer._flush()

    core.upsert_logs.assert_called_once()
    assert consumer._scalar_batch == []
    assert consumer._log_batch == []


def test_backoff_is_bounded(make_retry_consumer):
    consumer, core, clock = make_retry_consumer(retry_base_delay=1.0, retry_max_delay=2.0)
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("busy", retryable=True)

    delays = []
    for _ in range(6):
        consumer._flush()
        retry_at = consumer._retry_at
        assert retry_at is not None
        delays.append(retry_at - clock.now)
        clock.advance(retry_at - clock.now)

    assert delays[:3] == [1.0, 2.0, 2.0]
    assert max(delays) <= 2.0
    assert consumer._scalar_batch != []


@pytest.mark.parametrize(
    "retry_in, flush_timeout, expected",
    [
        (None, 0.5, 0.5),
        (10.0, 0.5, 0.5),
        (0.2, 0.5, 0.2),
        (0.0, 0.5, 0.0),
    ],
)
def test_queue_timeout_is_min_of_flush_timeout_and_remaining(make_retry_consumer, retry_in, flush_timeout, expected):
    consumer, _, clock = make_retry_consumer(flush_timeout=flush_timeout)
    consumer._retry_at = None if retry_in is None else clock.now + retry_in

    assert consumer._queue_timeout() == pytest.approx(expected)


def test_background_loop_retries_after_backoff(make_retry_consumer):
    """真实线程路径：失败 → queue.get 等待 → 自动重试成功。"""
    consumer, core, _ = make_retry_consumer(retry_base_delay=0.0, retry_max_delay=0.0, flush_timeout=0.05)
    record = _scalar()
    consumer._scalar_batch.append(record)
    failed = threading.Event()
    succeeded = threading.Event()

    def upsert(batch):
        if not failed.is_set():
            failed.set()
            raise CoreProtocolError("busy", retryable=True)
        succeeded.set()

    core.upsert_scalars.side_effect = upsert
    consumer.start()
    try:
        assert failed.wait(timeout=2.0), "consumer did not attempt the first flush"
        assert succeeded.wait(timeout=2.0), "consumer did not retry automatically"
    finally:
        consumer.stop()
        consumer._thread.join(timeout=2.0)

    assert not consumer._thread.is_alive()
    assert consumer._scalar_batch == []
    assert core.upsert_scalars.call_count == 2


def test_stop_with_pending_retry_keeps_records_and_warns(monkeypatch, make_retry_consumer):
    """停止时若仍处于退避窗口，保留记录并明确告警，不静默丢弃。"""
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer(retry_base_delay=5.0, flush_timeout=0.05)
    record = _scalar()
    consumer._scalar_batch.append(record)
    failed = threading.Event()

    def upsert(batch):
        failed.set()
        raise CoreProtocolError("busy", retryable=True)

    core.upsert_scalars.side_effect = upsert
    consumer.start()
    try:
        assert failed.wait(timeout=2.0), "consumer did not attempt the first flush"
        consumer.stop()
        consumer._thread.join(timeout=2.0)
        assert not consumer._thread.is_alive()
    finally:
        if consumer._thread.is_alive():
            consumer.stop()
            consumer._thread.join(timeout=2.0)

    assert consumer._scalar_batch == [record]
    assert len(warnings) == 1
    assert "did not accept 1 buffered record" in warnings[0][0]


def test_flush_chunks_batches_by_batch_size(make_retry_consumer):
    """单次 flush 按 batch_size 分片提交，避免积压后一次性发出超大请求。"""
    consumer, core, _ = make_retry_consumer(batch_size=2)
    records = [_scalar(key=f"loss/{i}") for i in range(5)]
    consumer._scalar_batch.extend(records)

    consumer._flush()

    assert [len(call.args[0]) for call in core.upsert_scalars.call_args_list] == [2, 2, 1]
    sent = [record for call in core.upsert_scalars.call_args_list for record in call.args[0]]
    assert sent == records
    assert consumer._scalar_batch == []


def test_retryable_failure_keeps_record_order_across_chunks(make_retry_consumer):
    """成功分片不回塞，失败分片与未发送尾巴按原顺序重试。"""
    consumer, core, clock = make_retry_consumer(batch_size=2)
    records = [_scalar(key=f"loss/{i}") for i in range(5)]
    consumer._scalar_batch.extend(records)
    core.upsert_scalars.side_effect = [None, CoreProtocolError("busy", retryable=True)]

    consumer._flush()

    assert consumer._scalar_batch == records[2:]
    assert [len(call.args[0]) for call in core.upsert_scalars.call_args_list] == [2, 2]

    core.upsert_scalars.side_effect = None
    core.upsert_scalars.reset_mock()
    clock.advance(0.5)
    consumer._flush()

    sent = [record for call in core.upsert_scalars.call_args_list for record in call.args[0]]
    assert sent == records[2:]
    assert consumer._scalar_batch == []


def test_terminal_failure_drops_only_rejected_chunk(monkeypatch, make_retry_consumer):
    """确定性拒绝只丢弃当前分片，剩余记录继续发送。"""
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer(batch_size=2)
    records = [_scalar(key=f"loss/{i}") for i in range(5)]
    consumer._scalar_batch.extend(records)
    core.upsert_scalars.side_effect = [None, CoreProtocolError("bad record", retryable=False), None]

    consumer._flush()

    assert consumer._scalar_batch == records[4:]
    assert len(warnings) == 1
    assert "count=2" in warnings[0][0]

    consumer._flush()

    sent = [record for call in core.upsert_scalars.call_args_list for record in call.args[0]]
    assert sent == records
    assert consumer._scalar_batch == []


def test_histogram_batch_is_restored_and_retried_on_retryable_failure(make_retry_consumer):
    """直方图走独立通道：可重试失败回塞 _histogram_batch，退避到期重提。"""
    consumer, core, clock = make_retry_consumer()
    record = _histogram()
    consumer._histogram_batch.append(record)
    core.upsert_histograms.side_effect = [CoreProtocolError("busy", retryable=True), None]

    consumer._flush()

    assert consumer._histogram_batch == [record]
    assert consumer._retry_due() is False

    clock.advance(0.5)
    consumer._flush()

    assert consumer._histogram_batch == []
    assert core.upsert_histograms.call_count == 2
    assert core.upsert_histograms.call_args_list[0].args[0] == [record]


def test_histogram_flush_notifies_callbacker(make_retry_consumer):
    consumer, core, _ = make_retry_consumer()
    callbacker = cast(MagicMock, consumer._callbacker)
    record = _histogram()
    consumer._histogram_batch.append(record)

    consumer._flush()

    core.upsert_histograms.assert_called_once_with([record])
    callbacker.on_histogram_flush.assert_called_once_with([record])


def test_terminal_histogram_error_drops_only_that_channel(monkeypatch, make_retry_consumer):
    """Go core 未实现 UpsertHistograms 时返回非可重试错误：只丢直方图，不影响其它通道。"""
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer()
    consumer._histogram_batch.append(_histogram())
    consumer._scalar_batch.append(_scalar())
    core.upsert_histograms.side_effect = CoreProtocolError("unimplemented", retryable=False)

    consumer._flush()

    assert consumer._histogram_batch == []
    assert consumer._scalar_batch == []
    core.upsert_scalars.assert_called_once()
    assert len(warnings) == 1
    assert "histograms" in warnings[0][0]


def test_columns_are_coalesced_before_chunking(make_retry_consumer):
    """列定义在分片前 coalesce，避免同一列跨分片重复上报。"""
    consumer, core, _ = make_retry_consumer(batch_size=2)
    first = ColumnRecord(column_key="a", column_type=ColumnType.COLUMN_TYPE_SCALAR, column_name="first")
    last = ColumnRecord(column_key="a", column_type=ColumnType.COLUMN_TYPE_SCALAR, column_name="last")
    other = ColumnRecord(column_key="b", column_type=ColumnType.COLUMN_TYPE_SCALAR, column_name="b")
    consumer._column_batch.extend([first, last, other])

    consumer._flush()

    assert [len(call.args[0]) for call in core.upsert_columns.call_args_list] == [2]
    sent = [record.column_name for call in core.upsert_columns.call_args_list for record in call.args[0]]
    assert sent == ["last", "b"]


def test_rpc_batch_size_chunks_scalars_larger_than_flush_threshold(make_retry_consumer):
    """标量按 rpc_batch_size 切片提交，logs 仍受较小 batch_size 约束。"""
    consumer, core, _ = make_retry_consumer(batch_size=2, rpc_batch_size=5)
    consumer._scalar_batch.extend([_scalar(key=f"s/{i}") for i in range(7)])
    consumer._log_batch.extend([object(), object(), object()])  # type: ignore[list-item]

    consumer._flush()

    assert [len(call.args[0]) for call in core.upsert_scalars.call_args_list] == [5, 2]
    assert [len(call.args[0]) for call in core.upsert_logs.call_args_list] == [2, 1]


def test_stop_drain_retries_while_progress_is_made(monkeypatch, make_retry_consumer):
    """退出排空期间 Core 短暂过载时持续重试直到全部刷完，且每次推进都会刷新停滞窗口。"""
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, clock = make_retry_consumer(
        batch_size=2,
        stop_timeout=5.0,
        retry_base_delay=1.0,
        retry_max_delay=2.0,
    )
    consumer._sleep = clock.advance
    records = [_scalar(key=f"loss/{i}") for i in range(5)]
    consumer._scalar_batch.extend(records)

    # 第 1 轮推进 2 条后过载；第 2 轮完全过载；第 3 轮推进 2 条后过载；第 4 轮完成剩余 1 条
    core.upsert_scalars.side_effect = [
        None,
        CoreProtocolError("busy", retryable=True),
        CoreProtocolError("busy", retryable=True),
        None,
        CoreProtocolError("busy", retryable=True),
        None,
    ]

    consumer._drain_before_stop()

    assert consumer._scalar_batch == []
    assert warnings == []
    sent = [rec for call in core.upsert_scalars.call_args_list for rec in call.args[0]]
    # 只统计成功接收的分片（索引 0, 3, 5）
    accepted = [rec for i in (0, 3, 5) for rec in core.upsert_scalars.call_args_list[i].args[0]]
    assert accepted == records
    assert len(sent) > len(records)


def test_stop_drain_gives_up_after_stall_timeout(make_retry_consumer):
    """退出排空期间 Core 持续零推进超过 stop_timeout 后停止重试并保留未接收记录。"""
    consumer, core, clock = make_retry_consumer(
        batch_size=2,
        stop_timeout=3.0,
        retry_base_delay=1.0,
        retry_max_delay=2.0,
    )
    consumer._sleep = clock.advance
    record = _scalar()
    consumer._scalar_batch.append(record)
    core.upsert_scalars.side_effect = CoreProtocolError("busy", retryable=True)

    consumer._drain_before_stop()

    assert consumer._scalar_batch == [record]
    assert clock.now == pytest.approx(103.0)
    assert core.upsert_scalars.call_count == 3


def test_stop_drain_waits_indefinitely_when_stop_timeout_is_none(monkeypatch, make_retry_consumer):
    """stop_timeout=None 时，只要 Core 返回可重试错误就持续阻塞重试，直到全部指标被接收。"""
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, clock = make_retry_consumer(
        batch_size=2,
        stop_timeout=None,
        retry_base_delay=1.0,
        retry_max_delay=5.0,
    )
    consumer._sleep = clock.advance
    record = _scalar()
    consumer._scalar_batch.append(record)
    # 连续失败 20 次（总等待时长远超 30s），第 21 次成功接收
    core.upsert_scalars.side_effect = [CoreProtocolError("busy", retryable=True)] * 20 + [None]

    consumer._drain_before_stop()

    assert consumer._scalar_batch == []
    assert warnings == []
    assert core.upsert_scalars.call_count == 21
    assert clock.now > 180.0


def test_stop_in_gather_branch_drains_until_all_records_accepted(make_retry_consumer):
    """生产默认 stop_timeout=None：STOP 落在聚拢分支时，持续重试直到积压全部提交。"""
    consumer, core, _ = make_retry_consumer(
        retry_base_delay=0.0,
        retry_max_delay=0.0,
        flush_timeout=0.05,
        stop_timeout=None,
    )
    for step in range(20):
        consumer._queue.put(MetricLogEvent(data={"loss": float(step)}, step=step, timestamp=Timestamp(seconds=step)))
    consumer.stop()
    core.upsert_scalars.side_effect = [
        CoreProtocolError("busy", retryable=True),
        CoreProtocolError("busy", retryable=True),
        None,
    ]

    consumer.start()
    try:
        consumer._thread.join(timeout=2.0)
    finally:
        if consumer._thread.is_alive():
            consumer._thread.join(timeout=2.0)

    assert not consumer._thread.is_alive()
    assert consumer._scalar_batch == []
    assert core.upsert_scalars.call_count == 3
    assert len(core.upsert_scalars.call_args_list[-1].args[0]) == 20
    callbacker = cast(MagicMock, consumer._callbacker)
    callbacker.on_scalar_flush.assert_called_once()
    assert len(callbacker.on_scalar_flush.call_args.args[0]) == 20


def test_stop_does_not_block_when_queue_full_and_consumer_exited(make_retry_consumer):
    full_queue: queue.Queue = queue.Queue(maxsize=1)
    full_queue.put(MetricLogEvent(data={"loss": 1.0}, step=1, timestamp=Timestamp()))
    consumer, _, _ = make_retry_consumer(event_queue=full_queue, stop_timeout=None)
    consumer._started = True  # noqa: SLF001

    consumer.stop()

    assert full_queue.qsize() == 1


def test_drain_before_stop_flushes_remaining_chunks_after_terminal_chunk(monkeypatch, make_retry_consumer):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    consumer, core, _ = make_retry_consumer(batch_size=2, stop_timeout=None)
    r1, r2, r3 = _scalar("loss/1"), _scalar("loss/2"), _scalar("loss/3")
    consumer._append_data_record(r1)
    consumer._append_data_record(r2)
    consumer._append_data_record(r3)
    core.upsert_scalars.side_effect = [CoreProtocolError("too large", retryable=False), None]

    consumer._drain_before_stop()

    assert consumer._batch_empty
    assert core.upsert_scalars.call_count == 2
    assert core.upsert_scalars.call_args_list[0].args[0] == [r1, r2]
    assert core.upsert_scalars.call_args_list[1].args[0] == [r3]
    assert len(warnings) == 1


def test_terminal_failure_warning_is_not_captured_as_user_log(monkeypatch, make_retry_consumer):
    """确定性拒绝的告警不能被终端代理采集，否则 core 拒收日志 → 告警 → 再上报形成自循环。"""
    seen = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: seen.append(pkg_console.capture_skipped()))
    consumer, core, _ = make_retry_consumer()
    consumer._scalar_batch.append(_scalar())
    core.upsert_scalars.side_effect = CoreProtocolError("bad record", retryable=False)

    consumer._flush()

    assert seen == [True]
    assert pkg_console.capture_skipped() is False


def test_buffer_cap_warning_is_not_captured_as_user_log(monkeypatch, make_retry_consumer):
    """缓冲超限告警同样不能被采集，否则丢弃告警本身又会变成一条待上传日志。"""
    seen = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: seen.append(pkg_console.capture_skipped()))
    consumer, _, _ = make_retry_consumer(max_buffered_records=1)

    consumer._append_data_record(_scalar(key="loss/0"))
    consumer._append_data_record(_scalar(key="loss/1"))

    assert seen == [True]
    assert pkg_console.capture_skipped() is False


def test_terminal_capture_log_rejection_does_not_self_loop(monkeypatch, make_retry_consumer):
    """闭环验证：StreamCapture → RunEmitter → BackgroundConsumer 在 core 持续拒绝日志时不自循环。"""
    stream = io.StringIO()
    monkeypatch.setattr(sys, "stdout", stream)
    emitter = RunEmitter()
    consumer, core, _ = make_retry_consumer()
    consumer._queue = emitter.queue
    core.upsert_logs.side_effect = CoreProtocolError("rejected log", retryable=False)

    def on_write(text: str) -> None:
        ts = Timestamp()
        ts.GetCurrentTime()
        emitter.emit(LogEvent(line=text, level=LogLevel.LOG_LEVEL_INFO, timestamp=ts))

    capture = StreamCapture("stdout", on_write, os.getpid())
    capture.install()
    try:
        print("user training log")
        rounds = 0
        while not emitter.queue.empty() and rounds < 5:
            rounds += 1
            while not emitter.queue.empty():
                consumer._handle_event(emitter.queue.get_nowait())
            consumer._flush()
    finally:
        capture.uninstall()

    screen = " ".join(stream.getvalue().split())
    assert rounds == 1
    assert emitter.queue.empty()
    assert core.upsert_logs.call_count == 1
    assert screen.count("Dropping logs records after terminal failure") == 1


def _too_large_error() -> CoreGoError:
    return CoreGoError(Stage.RPC, ReasonCode.REQUEST_TOO_LARGE, "request size 5242880 exceeds gRPC limit 4194304")


def test_oversized_chunk_is_split_until_only_the_offending_record_is_dropped(monkeypatch, make_retry_consumer):
    """单片超限按对半拆开重发：一条超限记录只丢自己，同片其它记录照常送达。"""
    consumer, core, _ = make_retry_consumer(batch_size=4)
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args[0]))
    sent = []

    def upsert_media(chunk):
        if any(record.key == "bad" for record in chunk):
            raise _too_large_error()
        sent.append([record.key for record in chunk])

    core.upsert_media.side_effect = upsert_media
    for key in ["a", "b", "bad", "c", "d", "e", "f", "g"]:
        consumer._media_batch.append(MediaRecord(key=key, step=1))

    consumer._flush()

    assert [key for chunk in sent for key in chunk] == ["a", "b", "c", "d", "e", "f", "g"]
    assert consumer._media_batch == []
    # 超限是确定性失败，不该触发退避重试
    assert consumer._retry_at is None
    assert len(warnings) == 1
    assert "'bad'" in warnings[0]


def test_single_oversized_record_is_dropped_with_its_key(monkeypatch, make_retry_consumer):
    """分片大小为 1 时无处可拆：只丢弃该条，并带指标 key 提示用户。"""
    consumer, core, _ = make_retry_consumer()
    warnings = []
    monkeypatch.setattr(
        pkg_console, "warning", lambda *args, **kwargs: warnings.append((args[0], pkg_console.capture_skipped()))
    )
    core.upsert_media.side_effect = _too_large_error()
    consumer._media_batch.append(MediaRecord(key="huge-image", step=3))

    consumer._flush()

    assert consumer._media_batch == []
    text, skipped = warnings[0]
    assert "huge-image" in text
    assert skipped is True


def test_normal_chunk_is_not_split(make_retry_consumer):
    """正常分片不触发二分：9 条记录按 batch_size=4 仍是 3 次 RPC。"""
    consumer, core, _ = make_retry_consumer(batch_size=4)
    for index in range(9):
        consumer._media_batch.append(MediaRecord(key=f"k{index}", step=1))

    consumer._flush()

    assert core.upsert_media.call_count == 3
