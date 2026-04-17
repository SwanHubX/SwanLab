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


def test_transport_put_warns_after_finish(make_scalar_record):
    """finish() 后 put() 静默丢弃并 warning。"""
    t = Transport(auto_start=False)
    t.start()
    t.finish()
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning") as mock_warning:
        t.put([make_scalar_record(step=1)])
    mock_warning.assert_called_once()


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
            return True, []

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

    t._thread.join.assert_called_once_with(timeout=Transport.FINISH_JOIN_TIMEOUT)
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
            return True, []

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
            return True, []

    t = Transport(auto_start=True)
    t._dispatcher = FakeDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1


# ─────────────────── consecutive failure ───────────────────


def test_transport_keeps_pending_records_and_warns_after_retry_exhaustion(make_scalar_record):
    """退避重试耗尽后应告警，但不丢 pending，新记录也保持在其后。"""
    attempts = []
    success_event = threading.Event()

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append([record.num for record in records])
            if len(attempts) == 1:
                return False, records
            success_event.set()
            return True, []

    t = Transport(batch_interval=0.01, auto_start=False)
    t._dispatcher = FlakyDispatch()  # type: ignore
    t.start()

    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 1
    second.num = 2

    with patch("swanlab.sdk.internal.core_python.transport.thread.console.warning") as mock_warning:
        t.put([first])
        time.sleep(0.05)
        t.put([second])
        assert success_event.wait(timeout=2.0)
        t.finish()

    assert attempts[0] == [1]
    assert attempts[1] == [1]
    assert attempts[2] == [2]
    mock_warning.assert_called_once()
