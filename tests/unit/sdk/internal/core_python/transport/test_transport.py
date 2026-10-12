import threading
import time
from unittest.mock import MagicMock, patch

import pytest

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.transport.worker import Transport
from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy, RetrySchedule
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult


def _ok_result(records: list[Record]) -> UploadResult:
    """构造确认成功的 UploadResult。"""
    return UploadResult(confirmed_records=list(records))


def _failed_result(records: list[Record], message: str = "upload failed") -> UploadResult:
    """构造失败保留的 UploadResult。"""
    return UploadResult(
        remaining_records=list(records),
        failures=[UploadFailure(source="swanlab_api", message=message, records=list(records))],
    )


def _fast_retry(*, base_delay: float = 0.01, max_delay: float = 0.04) -> RetryPolicy:
    """毫秒级退避策略，用于加速测试。"""
    return RetryPolicy(
        base_delay=base_delay,
        max_delay=max_delay,
        rng=lambda _low, _high: 1.0,
    )


class _RecordingPolicy(RetryPolicy):
    """记录 schedule 调用参数的重试策略替身。"""

    def __init__(self, *, base_delay: float = 0.01, max_delay: float = 0.04) -> None:
        super().__init__(
            base_delay=base_delay,
            max_delay=max_delay,
            rng=lambda _low, _high: 1.0,
        )
        self.calls: list[tuple[int, RetrySchedule]] = []

    def schedule(self, failures, failed_rounds: int) -> RetrySchedule:
        result = super().schedule(failures, failed_rounds)
        self.calls.append((failed_rounds, result))
        return result


def _wait_attempts(attempts: list, count: int, timeout: float = 2.0) -> bool:
    """轮询等待 dispatch 调用次数达到 count。"""
    deadline = time.monotonic() + timeout
    while len(attempts) < count and time.monotonic() < deadline:
        time.sleep(0.005)
    return len(attempts) >= count


def _make_transport(ctx, **kwargs) -> Transport:
    kwargs.setdefault("auto_start", False)
    return Transport(ctx=ctx, **kwargs)


def test_transport_start_creates_daemon_thread(mock_ctx):
    """start() 创建名为 THREAD_NAME 的守护线程。"""
    t = _make_transport(mock_ctx)
    t.start()
    try:
        assert t._thread is not None
        assert t._thread.daemon is True
        assert t._thread.name == Transport.THREAD_NAME
    finally:
        t.finish()


def test_transport_start_is_idempotent(mock_ctx):
    """重复调用 start() 不创建新线程。"""
    t = _make_transport(mock_ctx)
    t.start()
    first_thread = t._thread
    t.start()
    assert t._thread is first_thread
    t.finish()


def test_transport_start_uses_injected_sender_and_tracker(mock_ctx):
    """start() 将注入的 sender 与 tracker 绑定至 Dispatch。"""
    sender = MagicMock()
    tracker = MagicMock()
    t = _make_transport(mock_ctx, sender=sender, tracker=tracker)

    t.start()
    try:
        sender.set_tracker.assert_called_once_with(tracker)
        assert t._dispatcher is not None
        assert t._dispatcher._sender is sender
    finally:
        t.finish()


def test_transport_start_after_finish_is_noop(mock_ctx):
    """finish() 之后调用 start() 不创建新线程。"""
    t = _make_transport(mock_ctx)
    t.start()
    original_thread = t._thread
    t.finish()
    t.start()
    assert t._thread is original_thread


def test_transport_config_replacement_keeps_record_total(mock_ctx, make_config_record):
    tracker = UploadTracker()
    sender = MagicMock()
    sender.upload.side_effect = lambda _, records: _ok_result(records)
    t = _make_transport(mock_ctx, tracker=tracker, sender=sender)
    first = make_config_record()
    latest = make_config_record()
    first.num = latest.num = -3
    first.save.payload = b"first"
    latest.save.payload = b"latest"

    t.put([first])
    t.put([latest])
    assert tracker.snapshot().total_records == 1

    t.start()
    assert t.finish()
    sender.upload.assert_called_once_with("save", [latest])
    assert tracker.snapshot().uploaded_records == 1


def test_transport_uploads_latest_config_after_pending_retry(mock_ctx, make_config_record):
    tracker = UploadTracker()
    sender = MagicMock()
    t = _make_transport(mock_ctx, tracker=tracker, sender=sender, retry=_fast_retry())
    old = make_config_record()
    latest = make_config_record()
    old.num = latest.num = -3
    old.save.payload = b"old"
    latest.save.payload = b"latest"
    attempts = []
    latest_queued = threading.Event()

    def upload(_, records):
        attempts.append([record.save.payload for record in records])
        if len(attempts) == 1:
            t.put([latest])
            latest_queued.set()
            return _failed_result(records)
        return _ok_result(records)

    sender.upload.side_effect = upload
    t.put([old])
    t.start()
    try:
        assert latest_queued.wait(timeout=2.0)
        assert t.finish()
    finally:
        t.finish()

    assert attempts == [[b"old"], [b"old"], [b"latest"]]
    stats = tracker.snapshot()
    assert stats.total_records == stats.uploaded_records == 2


def test_transport_put_warns_after_finish(mock_ctx, make_scalar_record):
    """finish() 之后调用 put() 返回 False 并报错，空输入返回 True。"""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.error") as mock_error:
        assert t.put([make_scalar_record(step=1)]) is False
    mock_error.assert_called_once()
    accepted_transport = _make_transport(mock_ctx)
    assert accepted_transport.put([make_scalar_record(step=1)]) is True
    assert accepted_transport.put([]) is True


def test_transport_finish_is_idempotent(mock_ctx):
    """重复调用 finish() 保持幂等。"""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    t.finish()


def test_transport_finish_drains_remaining(make_ctx, make_scalar_record):
    """finish() 排空缓冲区剩余记录。"""
    dispatched = []

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            return _ok_result(records)

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx)
    t.start()
    t._dispatcher = FakeDispatch()  # type: ignore
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 1
    records[1].num = 2
    t.put(records)
    t.finish()
    assert len(dispatched) == 2


@pytest.mark.parametrize("alive", [False, True], ids=["drained", "still-running"])
def test_transport_finish_and_join_report_worker_state(mock_ctx, alive):
    """request_finish 不阻塞等待，finish/join 根据线程状态返回结果。"""
    t = _make_transport(mock_ctx)
    assert t.is_alive() is False
    t._thread = MagicMock()
    t._thread.is_alive.return_value = alive
    t._drained = not alive
    assert t.is_alive() is alive

    t.request_finish()

    assert t._finished is True
    t._thread.join.assert_not_called()

    assert t.finish() is (not alive)
    t._thread.join.assert_called_once_with(timeout=Transport.FINISH_JOIN_TIMEOUT)
    t._thread.join.reset_mock()
    assert t.join(timeout=0) is (not alive)
    t._thread.join.assert_called_once_with(timeout=0)


def test_transport_marks_tracker_finished_when_worker_exits(mock_ctx):
    tracker = UploadTracker()
    tracker.set_state(CoreState.CORE_STATE_RUNNING)
    t = _make_transport(mock_ctx, tracker=tracker)

    t.start()
    t.request_finish()

    assert t.join(timeout=2.0) is True
    assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED


@pytest.mark.parametrize("auto_start", [False, True], ids=["manual_start", "auto_start"])
def test_transport_drains_on_batch_interval(make_ctx, make_scalar_record, auto_start):
    """put() 后到达 batch_interval 时排空并分发记录。"""
    dispatched = []
    event = threading.Event()

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            event.set()
            return _ok_result(records)

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, auto_start=auto_start)
    if not auto_start:
        t.start()
    t._dispatcher = FakeDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1


def test_transport_outage_recovers_and_resumes_backlog(make_ctx, make_scalar_record):
    """服务停机期间接受新记录入队并保留未确认记录，恢复后完成续传与排空。"""
    sender = MagicMock()
    tracker = UploadTracker()
    attempts: list[tuple[list[int], str]] = []
    first_failure_seen = threading.Event()
    recovered = threading.Event()

    def upload(_record_type, records):
        nums = [record.num for record in records]
        if not recovered.is_set():
            attempts.append((nums, "fail"))
            first_failure_seen.set()
            return _failed_result(records, message="http 503")
        attempts.append((nums, "ok"))
        return _ok_result(records)

    sender.upload.side_effect = upload
    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, sender=sender, tracker=tracker, retry=_fast_retry())
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num, second.num = 1, 2

    t.start()
    with patch("swanlab.sdk.internal.core_python.transport.worker.UploadWarningThrottle.warn") as mock_warn:
        t.put([first])
        try:
            assert first_failure_seen.wait(timeout=2.0)
            deadline = time.monotonic() + 2.0
            while not t.last_failures() and time.monotonic() < deadline:
                time.sleep(0.005)
            failure = t.last_failures()[0]
            assert failure.message == "http 503"
            assert failure.records == [first]

            assert t.put([second]) is True
            assert t.is_alive() is True
            pending_nums = sorted(record.num for record in t.pending_records())
            assert pending_nums == [1, 2]
            assert tracker.snapshot().total_records == 2

            recovered.set()
            assert t.finish()

            failed_attempts = [attempt for attempt in attempts if attempt[1] == "fail"]
            assert failed_attempts
            assert all(attempt == ([1], "fail") for attempt in failed_attempts)
            ok_attempts = [attempt for attempt in attempts if attempt[1] == "ok"]
            assert ok_attempts == [([1], "ok"), ([2], "ok")]
            assert mock_warn.called
            stats = tracker.snapshot()
            assert stats.uploaded_records == stats.total_records == 2
            assert t.pending_records() == []
            assert t.last_failures() == []
            sender.close.assert_called_once()
        finally:
            recovered.set()
            t.finish()


def test_transport_retry_waits_for_deadline_and_wake_does_not_bypass(make_ctx, make_scalar_record):
    """失败重试等待退避截止时间，request_finish 唤醒不跳过等待。"""
    attempts = []
    first_dispatch_done = threading.Event()

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append(time.monotonic())
            if len(attempts) == 1:
                first_dispatch_done.set()
                return _failed_result(records)
            return _ok_result(records)

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, retry=RetryPolicy(base_delay=0.5, max_delay=0.5, rng=lambda _low, _high: 1.0))
    t.start()
    t._dispatcher = FlakyDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert first_dispatch_done.wait(timeout=2.0)

    t.request_finish()
    time.sleep(0.15)
    assert len(attempts) == 1

    assert t.join(timeout=3.0), "worker should retry after the backoff deadline and finish"
    assert len(attempts) == 2
    assert attempts[1] - attempts[0] >= 0.4


def test_transport_backoff_counter_resets_after_success_round(make_ctx, make_scalar_record):
    """连续失败轮数在整轮成功后归零。"""
    attempts = []
    second = make_scalar_record(step=2)
    second.num = 2

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append(len(attempts))
            if len(attempts) == 3:
                t.put([second])
                return _ok_result(records)
            if len(attempts) <= 2 or len(attempts) == 4:
                return _failed_result(records)
            return _ok_result(records)

    policy = _RecordingPolicy()
    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, retry=policy)
    t.start()
    t._dispatcher = FlakyDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert _wait_attempts(attempts, 4)
    assert t.finish()
    assert len(attempts) == 5

    assert [round_ for round_, _ in policy.calls] == [1, 2, 1]
    delays = [schedule.delay for _, schedule in policy.calls]
    assert delays[0] < delays[1]
    assert delays[2] == delays[0]


def test_transport_blocked_failure_uses_blocked_interval(make_ctx, make_scalar_record):
    """local_data 失败使用 blocked 间隔调度，恢复后完成排空。"""
    attempts = []

    class BlockedDispatch:
        def __call__(self, records):
            attempts.append(len(attempts))
            if len(attempts) == 1:
                return UploadResult(
                    remaining_records=list(records),
                    failures=[UploadFailure(source="local_data", message="save file not found", records=list(records))],
                )
            return _ok_result(records)

    policy = _RecordingPolicy(base_delay=0.05, max_delay=0.05)
    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, retry=policy)
    t.start()
    t._dispatcher = BlockedDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert _wait_attempts(attempts, 2)
    assert t.finish()

    assert len(policy.calls) == 1
    assert policy.calls[0][1].delay == 0.05
    assert t.pending_records() == []


def test_transport_huge_retry_delay_waits_in_slices_without_crashing(make_ctx, make_scalar_record, monkeypatch):
    """大重试间隔按 TIMEOUT_MAX 分段等待。"""
    attempts = []
    waiting = threading.Event()
    timeout_max = 0.05
    monkeypatch.setattr(threading, "TIMEOUT_MAX", timeout_max)

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append(len(attempts))
            if len(attempts) == 1:
                return _failed_result(records)
            return _ok_result(records)

    class HugeDelayPolicy(RetryPolicy):
        def schedule(self, failures, failed_rounds: int) -> RetrySchedule:
            return RetrySchedule(1e20)

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, retry=HugeDelayPolicy())
    original_wait = t._cond.wait

    def wait_in_slice(timeout=None):
        assert timeout is not None and 0 <= timeout <= timeout_max
        if timeout == timeout_max:
            waiting.set()
        return original_wait(timeout=timeout)

    monkeypatch.setattr(t._cond, "wait", wait_in_slice)
    t.start()
    t._dispatcher = FlakyDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1

    t.put([record])
    try:
        assert waiting.wait(timeout=2.0)
        assert t.is_alive() is True
        assert len(attempts) == 1
        assert t.pending_records() == [record]

        with t._cond:
            t._retry_at = t._retry.clock()
            t._cond.notify_all()
        assert t.finish(timeout=2.0)
        assert len(attempts) == 2
    finally:
        with t._cond:
            t._retry_at = t._retry.clock()
            t._cond.notify_all()
        t.finish(timeout=2.0)


def test_transport_worker_crash_retains_data_and_finish_reports_failure(make_ctx, make_scalar_record):
    """worker 内部异常退出时保留未确认数据，finish() 返回 False。"""
    record = make_scalar_record(step=1)
    record.num = 1
    crashed = threading.Event()
    tracker = UploadTracker()
    sender = MagicMock()

    class CrashingDispatch:
        def __call__(self, records):
            crashed.set()
            return object()

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, tracker=tracker, sender=sender, retry=_fast_retry())
    t.start()
    t._dispatcher = CrashingDispatch()  # type: ignore
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.error") as mock_error:
        t.put([record])
        assert crashed.wait(timeout=2.0)
        assert t.join(timeout=3.0), "worker should exit after unexpected error"
        assert any("unexpectedly" in str(call.args[0]) for call in mock_error.call_args_list)

    assert t.finish() is False
    assert t.pending_records() == [record]
    assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED
    sender.close.assert_called_once()
    summary = t.finish_summary()
    assert "worker crashed: AttributeError" in summary
    assert "1 record(s) retained" in summary


def test_transport_finish_summarizes_last_failure_and_fallback(mock_ctx, make_scalar_record):
    """finish_summary() 包含保留记录数与最近失败摘要。"""
    t = _make_transport(mock_ctx)
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    with t._cond:
        t._last_failures = [UploadFailure(source="swanlab_api", message="http 503", status=503)]

    summary = t.finish_summary()

    assert "1 record(s) retained" in summary
    assert "last failure: swanlab_api - http 503" in summary

    with t._cond:
        t._last_failures = []
    assert "worker exited unexpectedly" in t.finish_summary()


def test_transport_closes_sender_on_clean_drain(mock_ctx):
    """排空退出时关闭 sender 资源。"""
    sender = MagicMock()
    t = _make_transport(mock_ctx, sender=sender)
    t.start()
    assert t.finish()
    sender.close.assert_called_once()


def test_transport_finish_false_when_never_started_with_buffered_records(mock_ctx, make_scalar_record):
    """未启动且缓冲区非空时 finish() 返回 False。"""
    t = _make_transport(mock_ctx)
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert t.finish() is False

    empty = _make_transport(mock_ctx)
    assert empty.finish() is True
