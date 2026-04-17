import threading
from unittest.mock import patch

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
    t.put(records)
    assert len(t._buffer) == 1


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
    t.put(records)
    t.finish()
    assert len(dispatched) == 2


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
    t.put([make_scalar_record(step=1)])
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
    t.put([make_scalar_record(step=1)])
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
                t._buffer.extend(records)
            if call_count >= 3:
                enough_failures.set()
            return False

    t = Transport(batch_interval=0.01, auto_start=False)
    t._dispatcher = FailingDispatch()  # type: ignore
    t.start()
    t.put([make_scalar_record(step=1)])
    enough_failures.wait(timeout=5)
    t.finish()
    assert t._thread is not None
    t._thread.join(timeout=5)
    assert not t._thread.is_alive()
    assert call_count >= 3
