import threading
import time
from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.core_python.transport.thread import Transport

# ─────────────────── 默认值 ───────────────────


def test_transport_defaults():
    """默认 batch_interval = 1.0。"""
    t = Transport(auto_start=False)
    assert t._batch_interval == Transport.BATCH_INTERVAL
    assert t._batch_interval == 1.0


def test_transport_custom_batch_interval():
    """自定义 batch_interval。"""
    t = Transport(batch_interval=0.5, auto_start=False)
    assert t._batch_interval == 0.5


# ─────────────────── start ───────────────────


def test_transport_start_creates_daemon_thread():
    """start() 创建守护线程，名称为 THREAD_NAME。"""
    t = Transport(auto_start=False)
    t.start()
    try:
        assert t._thread is not None
        assert t._thread.daemon is True
        assert t._thread.name == Transport.THREAD_NAME
    finally:
        t.finish()


def test_transport_start_is_idempotent():
    """重复 start() 不创建第二个线程。"""
    t = Transport(auto_start=False)
    t.start()
    first_thread = t._thread
    t.start()
    assert t._thread is first_thread
    t.finish()


def test_transport_start_after_finish_is_noop():
    """finish() 后 start() 无效。"""
    t = Transport(auto_start=False)
    t.start()
    t.finish()
    t.start()
    assert t._thread is not None
    # thread 已 join，不应创建新线程


# ─────────────────── put ───────────────────


def test_transport_put_appends_to_buffer(make_scalar_record):
    """put() 后 buffer 增长。"""
    t = Transport(auto_start=False)
    records = [make_scalar_record(step=1)]
    records[0].num = 1
    t.put(records)
    assert len(t._buf) == 1


def test_transport_put_dedups_records_by_num(make_scalar_record):
    """已编号 record 视为稳定事件，按 num 去重。"""
    t = Transport(auto_start=False)
    record = make_scalar_record(step=1)
    record.num = 7

    t.put([record])
    t.put([record])

    assert t._buf.drain() == [record]


def test_transport_put_keeps_distinct_record_nums(make_scalar_record):
    """不同 num 的 record 应同时保留。"""
    t = Transport(auto_start=False)
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=1)
    first.num = 8
    second.num = 9

    t.put([first, second])

    assert t._buf.drain() == [first, second]


def test_transport_put_raises_after_finish(make_scalar_record):
    """finish() 后 put() 抛 RuntimeError。"""
    t = Transport(auto_start=False)
    t.start()
    t.finish()
    with patch.object(t, "_finished", True):
        try:
            t.put([make_scalar_record(step=1)])
        except RuntimeError:
            pass


# ─────────────────── finish ───────────────────


def test_transport_finish_is_idempotent():
    """重复 finish() 只执行一次。"""
    t = Transport(auto_start=False)
    t.start()
    t.finish()
    t.finish()  # should not raise


def test_transport_finish_drains_remaining(make_scalar_record):
    """finish() 排空残余 buffer。"""
    dispatched = []

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            return True

    t = Transport(auto_start=False)
    t._dispatcher = FakeDispatch()  # type: ignore
    t.start()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 1
    records[1].num = 2
    t.put(records)
    t.finish()
    assert len(dispatched) == 2


def test_transport_finish_does_not_close_sender_before_thread_stops():
    """若 join 超时且线程仍存活，finish() 不应提前关闭 sender。"""
    sender = MagicMock()
    t = Transport(sender=sender, auto_start=False)
    t._thread = MagicMock()
    t._thread.is_alive.return_value = True

    t.finish()

    t._thread.join.assert_called_once_with(timeout=5)
    sender.close.assert_not_called()


# ─────────────────── 事件驱动 ───────────────────


def test_transport_put_wakes_thread_immediately(make_scalar_record):
    """put() 后线程立刻唤醒，不等 timeout。"""
    dispatched = []
    event = threading.Event()

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            event.set()
            return True

    t = Transport(batch_interval=10.0, auto_start=False)
    t._dispatcher = FakeDispatch()  # type: ignore
    t.start()
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    # batch_interval=10 但 put 应立刻唤醒，0.5s 足够验证
    event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1


# ─────────────────── auto_start 端到端 ───────────────────


def test_transport_integration_auto_start(make_scalar_record):
    """auto_start=True 时自动启动并正常处理。"""
    dispatched = []
    event = threading.Event()

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            event.set()
            return True

    t = Transport(auto_start=True)
    t._dispatcher = FakeDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1


# ─────────────────── consecutive failure ───────────────────


def test_transport_drops_after_max_consecutive_failures(make_scalar_record):
    """连续回滚超限后，finish 时放弃 buffer 退出。"""
    call_count = 0
    enough_failures = threading.Event()

    class FailingDispatch:
        def __call__(self, records):
            nonlocal call_count
            call_count += 1
            with t._cond:
                t._buf.prepend(records)
            if call_count >= 3:
                enough_failures.set()
            return False

    t = Transport(batch_interval=0.01, auto_start=False)
    t._dispatcher = FailingDispatch()  # type: ignore
    t.start()
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    enough_failures.wait(timeout=5)
    t.finish()
    assert t._thread is not None
    t._thread.join(timeout=5)
    assert not t._thread.is_alive()
    assert call_count >= 3


def test_transport_finish_interrupts_failure_backoff(make_scalar_record):
    """finish() 应打断失败退避，而不是被 sleep 卡住。"""
    first_failure = threading.Event()

    class FailingDispatch:
        def __call__(self, records):
            with t._cond:
                t._buf.prepend(records)
            first_failure.set()
            return False

    t = Transport(batch_interval=1.0, auto_start=False)
    t._dispatcher = FailingDispatch()  # type: ignore
    t.start()

    try:
        record = make_scalar_record(step=1)
        record.num = 1
        t.put([record])
        assert first_failure.wait(timeout=2.0)

        started_at = time.monotonic()
        t.finish()
        elapsed = time.monotonic() - started_at

        assert elapsed < 0.5
    finally:
        if not t._finished:
            t.finish()
