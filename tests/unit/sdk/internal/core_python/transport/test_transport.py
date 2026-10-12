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
    """构造全部确认的 Dispatch 结果（测试替身用）。"""
    return UploadResult(confirmed_records=list(records))


def _failed_result(records: list[Record], message: str = "upload failed") -> UploadResult:
    """构造全部保留的 Dispatch 失败结果（测试替身用）。"""
    return UploadResult(
        remaining_records=list(records),
        failures=[UploadFailure(source="swanlab_api", message=message, records=list(records))],
    )


def _fast_retry(*, base_delay: float = 0.01, max_delay: float = 0.04, blocked_delay: float = 0.02) -> RetryPolicy:
    """毫秒级退避的极速策略：避免测试真实等待，抖动固定 1.0 保证确定性。"""
    return RetryPolicy(
        base_delay=base_delay,
        max_delay=max_delay,
        blocked_delay=blocked_delay,
        rng=lambda _low, _high: 1.0,
    )


class _RecordingPolicy(RetryPolicy):
    """记录 schedule 调用的策略替身：捕获 (failed_rounds, schedule) 序列供断言。"""

    def __init__(self, *, base_delay: float = 0.01, max_delay: float = 0.04, blocked_delay: float = 0.02) -> None:
        super().__init__(
            base_delay=base_delay,
            max_delay=max_delay,
            blocked_delay=blocked_delay,
            rng=lambda _low, _high: 1.0,
        )
        self.calls: list[tuple[int, RetrySchedule]] = []

    def schedule(self, failures, failed_rounds: int) -> RetrySchedule:
        result = super().schedule(failures, failed_rounds)
        self.calls.append((failed_rounds, result))
        return result


def _wait_attempts(attempts: list, count: int, timeout: float = 2.0) -> bool:
    """轮询等待 dispatch 尝试次数到达 count（事件在 dispatch 内 set 早于结果应用，轮询更稳）。"""
    deadline = time.monotonic() + timeout
    while len(attempts) < count and time.monotonic() < deadline:
        time.sleep(0.005)
    return len(attempts) >= count


def _make_transport(ctx, **kwargs) -> Transport:
    kwargs.setdefault("auto_start", False)
    return Transport(ctx=ctx, **kwargs)


# ─────────────────── defaults ───────────────────


@pytest.mark.parametrize("custom_interval", [None, 0.5], ids=["default", "custom"])
def test_transport_batch_interval_reads_ctx_settings(make_ctx, custom_interval):
    """batch_interval 从 ctx 设置读取：默认与自定义值均生效。"""
    overrides = {} if custom_interval is None else {"batch_interval": custom_interval}
    ctx = make_ctx(**overrides)
    t = _make_transport(ctx)
    assert t._batch_interval == ctx.config.record_interval
    if custom_interval is not None:
        assert t._batch_interval == custom_interval


# ─────────────────── start ───────────────────


def test_transport_start_creates_daemon_thread(mock_ctx):
    """start() creates a daemon thread named THREAD_NAME."""
    t = _make_transport(mock_ctx)
    t.start()
    try:
        assert t._thread is not None
        assert t._thread.daemon is True
        assert t._thread.name == Transport.THREAD_NAME
    finally:
        t.finish()


def test_transport_start_is_idempotent(mock_ctx):
    """Repeated start() does not create a second thread."""
    t = _make_transport(mock_ctx)
    t.start()
    first_thread = t._thread
    t.start()
    assert t._thread is first_thread
    t.finish()


def test_transport_start_uses_injected_sender_and_tracker(mock_ctx):
    """start() wires an injected sender to Dispatch and injects the tracker into it."""
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
    """start() after finish() is a no-op：不创建新线程，线程身份保持不变。"""
    t = _make_transport(mock_ctx)
    t.start()
    original_thread = t._thread
    t.finish()
    t.start()
    assert t._thread is original_thread


# ─────────────────── put ───────────────────


def test_transport_put_dedups_records_by_num(mock_ctx, make_scalar_record):
    """Records with the same num are deduplicated."""
    t = _make_transport(mock_ctx)
    record = make_scalar_record(step=1)
    record.num = 7

    t.put([record])
    t.put([record])

    assert t._buf.drain() == [record]


def test_transport_put_keeps_distinct_record_nums(mock_ctx, make_scalar_record):
    """Records with different nums are kept separately."""
    t = _make_transport(mock_ctx)
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=1)
    first.num = 8
    second.num = 9

    t.put([first, second])

    assert t._buf.drain() == [first, second]


def test_transport_put_tracks_queued_record_total_by_default(mock_ctx, make_scalar_record):
    tracker = UploadTracker()
    t = _make_transport(mock_ctx, tracker=tracker)
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 1
    records[1].num = 2

    t.put(records)

    assert tracker.snapshot().total_records == 2


def test_transport_put_tracks_only_deduped_records(mock_ctx, make_scalar_record):
    tracker = UploadTracker()
    t = _make_transport(mock_ctx, tracker=tracker)
    record = make_scalar_record(step=1)
    record.num = 7

    t.put([record])
    t.put([record])

    assert tracker.snapshot().total_records == 1


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
    # 统计推进由 Dispatch 依据 confirmed_records 完成
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
    """put() after finish() 在锁内判定拒绝：返回 False 并报错，不只打印（plan §5）。"""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.error") as mock_error:
        assert t.put([make_scalar_record(step=1)]) is False
    mock_error.assert_called_once()
    # 接受路径返回 True；空输入视为无操作（非拒绝）
    accepted_transport = _make_transport(mock_ctx)
    assert accepted_transport.put([make_scalar_record(step=1)]) is True
    assert accepted_transport.put([]) is True


# ─────────────────── finish ───────────────────


def test_transport_finish_is_idempotent(mock_ctx):
    """Repeated finish() only runs once."""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    t.finish()  # should not raise


def test_transport_finish_drains_remaining(make_ctx, make_scalar_record):
    """finish() drains the remaining buffer."""
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


def test_transport_finish_waits_for_thread(mock_ctx):
    """finish() calls join with timeout from class constant; 排空成功才返回 True。"""
    t = _make_transport(mock_ctx)
    t._thread = MagicMock()
    t._thread.is_alive.return_value = False
    t._drained = True  # 模拟线程已正常排空退出

    assert t.finish() is True

    t._thread.join.assert_called_once_with(timeout=Transport.FINISH_JOIN_TIMEOUT)


def test_transport_request_finish_does_not_join_thread(mock_ctx):
    """request_finish() only notifies the worker and returns without joining."""
    t = _make_transport(mock_ctx)
    t._thread = MagicMock()

    t.request_finish()

    assert t._finished is True
    t._thread.join.assert_not_called()


def test_transport_marks_tracker_finished_when_worker_exits(mock_ctx):
    tracker = UploadTracker()
    tracker.set_state(CoreState.CORE_STATE_RUNNING)
    t = _make_transport(mock_ctx, tracker=tracker)

    t.start()
    t.request_finish()

    assert t.join(timeout=2.0) is True
    assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED


def test_transport_join_reports_timeout(mock_ctx):
    """join() returns False when the worker thread is still alive after timeout."""
    t = _make_transport(mock_ctx)
    t._thread = MagicMock()
    t._thread.is_alive.return_value = True

    assert t.join(timeout=0) is False
    t._thread.join.assert_called_once_with(timeout=0)


def test_transport_is_alive_reflects_thread_state(mock_ctx):
    t = _make_transport(mock_ctx)
    assert t.is_alive() is False

    t._thread = MagicMock()
    t._thread.is_alive.return_value = True

    assert t.is_alive() is True


# ─────────────────── 定时攒批 ───────────────────


@pytest.mark.parametrize("auto_start", [False, True], ids=["manual_start", "auto_start"])
def test_transport_drains_on_batch_interval(make_ctx, make_scalar_record, auto_start):
    """After put(), the thread auto-drains and dispatches on batch_interval timeout.
    Covers both manual start() and auto_start=True constructor paths.
    """
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
    """停机恢复：503 持续失败期间可继续入队、记录与失败原因可查询、节流告警触发、worker 不退出；
    恢复后旧 pending 与停机期间的新记录全部续传，确认计数正确。

    恢复时机由测试线程经 recovered 事件控制：设置前 sender 恒返回失败，
    停机期间的断言不依赖测试线程的调度速度；失败尝试次数随调度变化，
    只断言恢复后旧记录、新记录各成功一次。
    """
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
            # 失败事件在 dispatch 内 set，早于结果应用到实例状态，这里轮询等待登记完成
            deadline = time.monotonic() + 2.0
            while not t.last_failures() and time.monotonic() < deadline:
                time.sleep(0.005)
            failure = t.last_failures()[0]
            assert failure.message == "http 503"
            assert failure.records == [first]

            # 停机期间：新记录入队成功，worker 存活，两条记录均可查询
            assert t.put([second]) is True
            assert t.is_alive() is True
            pending_nums = sorted(record.num for record in t.pending_records())
            assert pending_nums == [1, 2]
            assert tracker.snapshot().total_records == 2

            recovered.set()
            assert t.finish()

            # 停机期间仅旧记录反复失败（次数随调度变化，不小于一次）；
            # 恢复后旧 pending 先续传成功，随后新记录排空，各恰好成功一次
            failed_attempts = [attempt for attempt in attempts if attempt[1] == "fail"]
            assert failed_attempts
            assert all(attempt == ([1], "fail") for attempt in failed_attempts)
            ok_attempts = [attempt for attempt in attempts if attempt[1] == "ok"]
            assert ok_attempts == [([1], "ok"), ([2], "ok")]
            # 失败轮触发节流告警（每个失败轮调用一次 warn，实际打印频率由节流器控制）
            assert mock_warn.called
            stats = tracker.snapshot()
            assert stats.uploaded_records == stats.total_records == 2
            assert t.pending_records() == []
            assert t.last_failures() == []
            sender.close.assert_called_once()
        finally:
            # 断言失败也先解除停机再关闭 worker，避免遗留自旋线程
            recovered.set()
            t.finish()


# ─────────────────── 退避重试调度（P3） ───────────────────


def test_transport_retry_waits_for_deadline_and_wake_does_not_bypass(make_ctx, make_scalar_record):
    """失败后按退避截止时间重试：唤醒（finish 通知）不绕过截止时间，不提前重试。"""
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
    # 真实时钟 + 0.5s 首次退避：足够区分"立即重试"与"按截止时间重试"
    t = _make_transport(ctx, retry=RetryPolicy(base_delay=0.5, max_delay=0.5, rng=lambda _low, _high: 1.0))
    t.start()
    t._dispatcher = FlakyDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert first_dispatch_done.wait(timeout=2.0)

    t.request_finish()  # 唤醒 worker：必须重新检查截止时间，不得提前重试
    time.sleep(0.15)
    assert len(attempts) == 1

    assert t.join(timeout=3.0), "worker should retry after the backoff deadline and finish"
    assert len(attempts) == 2
    # 首次重试发生在退避截止时间之后（0.5s ± 调度误差）
    assert attempts[1] - attempts[0] >= 0.4


def test_transport_backoff_counter_resets_after_success_round(make_ctx, make_scalar_record):
    """连续失败轮计数递增（1、2），整轮成功后清零，再次失败从 1 重新开始（plan §6.2）。"""
    attempts = []
    second = make_scalar_record(step=2)
    second.num = 2

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append(len(attempts))
            if len(attempts) == 3:
                # 本轮成功，同时补一条新记录驱动第 4 轮
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
    # 失败(1) → 失败(2) → 成功(重置) → 失败(1) → 成功排空
    assert _wait_attempts(attempts, 4)
    assert t.finish()
    assert len(attempts) == 5

    assert [round_ for round_, _ in policy.calls] == [1, 2, 1]
    # 退避随计数指数增长后回落
    delays = [schedule.delay for _, schedule in policy.calls]
    assert delays[0] < delays[1]
    assert delays[2] == delays[0]


def test_transport_blocked_failure_uses_blocked_interval(make_ctx, make_scalar_record):
    """本地坏数据等确定性失败走 blocked 低频间隔，恢复后仍可排空（文件恢复可自愈）。"""
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

    policy = _RecordingPolicy(blocked_delay=0.05)
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
    assert policy.calls[0][1].kind == "blocked"
    assert policy.calls[0][1].delay == 0.05
    assert t.pending_records() == []


def test_transport_huge_retry_delay_waits_in_slices_without_crashing(make_ctx, make_scalar_record, monkeypatch):
    """超大重试间隔按平台上限分段等待，保留记录且不提前重试。"""
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
            return RetrySchedule("rate_limited", 1e20)

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


# ─────────────────── finish 语义与资源关闭（P4） ───────────────────


def test_transport_worker_crash_retains_data_and_finish_reports_failure(make_ctx, make_scalar_record):
    """worker 内部未预期异常：线程退出、sender 关闭、tracker 置 FINISHED，
    但 finish 为 False 且数据保留——不伪报成功（plan §8.4）。"""
    record = make_scalar_record(step=1)
    record.num = 1
    crashed = threading.Event()
    tracker = UploadTracker()
    sender = MagicMock()

    class CrashingDispatch:
        def __call__(self, records):
            crashed.set()
            return object()  # 非 UploadResult：worker 在结果解包处崩溃

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, tracker=tracker, sender=sender, retry=_fast_retry())
    t.start()
    t._dispatcher = CrashingDispatch()  # type: ignore
    with patch("swanlab.sdk.internal.core_python.transport.worker.console.error") as mock_error:
        t.put([record])
        assert crashed.wait(timeout=2.0)
        assert t.join(timeout=3.0), "worker should exit after unexpected error"
        assert any("unexpectedly" in str(call.args[0]) for call in mock_error.call_args_list)

    assert t.finish() is False  # 线程已退出但未排空，不得伪报成功
    assert t.pending_records() == [record]
    assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED  # 进度轮询可退出
    sender.close.assert_called_once()  # 崩溃路径同样关闭 sender 资源
    # 摘要优先报告崩溃类型并携带积压数量，不拼接原始异常文本
    summary = t.finish_summary()
    assert "worker crashed: AttributeError" in summary
    assert "1 record(s) retained" in summary


def test_transport_finish_summarizes_last_failure_and_fallback(mock_ctx, make_scalar_record):
    """无崩溃时：最近一轮失败进入摘要，无失败回退通用文案；积压数量含 buffer 与 pending。"""
    t = _make_transport(mock_ctx)  # 从未启动
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
    """正常排空退出时关闭 sender 自持资源。"""
    sender = MagicMock()
    t = _make_transport(mock_ctx, sender=sender)
    t.start()
    assert t.finish()
    sender.close.assert_called_once()


def test_transport_finish_false_when_never_started_with_buffered_records(mock_ctx, make_scalar_record):
    """从未启动但缓冲有数据：finish 不得报告成功。"""
    t = _make_transport(mock_ctx)  # auto_start=False
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert t.finish() is False

    empty = _make_transport(mock_ctx)
    assert empty.finish() is True
