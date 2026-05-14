import threading
from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.core_python.transport.thread import Transport


def _make_transport(ctx, **kwargs) -> Transport:
    kwargs.setdefault("sender", MagicMock())
    kwargs.setdefault("auto_start", False)
    return Transport(ctx=ctx, **kwargs)


# ─────────────────── defaults ───────────────────


def test_transport_defaults(mock_ctx):
    """Default batch_interval reads from ctx settings."""
    t = _make_transport(mock_ctx)
    assert t._batch_interval == mock_ctx.config.record_interval


def test_transport_custom_batch_interval(make_ctx):
    """Custom batch_interval via ctx settings."""
    ctx = make_ctx(batch_interval=0.5)
    t = _make_transport(ctx)
    assert t._batch_interval == 0.5


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


def test_transport_start_after_finish_is_noop(mock_ctx):
    """start() after finish() is a no-op."""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    t.start()
    assert t._thread is not None


# ─────────────────── put ───────────────────


def test_transport_put_appends_to_buffer(mock_ctx, make_scalar_record):
    """put() grows the buffer."""
    t = _make_transport(mock_ctx)
    records = [make_scalar_record(step=1)]
    records[0].num = 1
    t.put(records)
    assert len(t._buf) == 1


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


def test_transport_put_warns_after_finish(mock_ctx, make_scalar_record):
    """put() after finish() logs an error."""
    t = _make_transport(mock_ctx)
    t.start()
    t.finish()
    with patch("swanlab.sdk.internal.core_python.transport.thread.console.error") as mock_error:
        t.put([make_scalar_record(step=1)])
    mock_error.assert_called_once()


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
            return True, []

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx)
    t._dispatcher = FakeDispatch()  # type: ignore
    t.start()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    records[0].num = 1
    records[1].num = 2
    t.put(records)
    t.finish()
    assert len(dispatched) == 2


def test_transport_finish_waits_for_thread(mock_ctx):
    """finish() calls join with timeout from class constant."""
    sender = MagicMock()
    t = _make_transport(mock_ctx, sender=sender)
    t._thread = MagicMock()

    t.finish()

    t._thread.join.assert_called_once_with(timeout=Transport.FINISH_JOIN_TIMEOUT)


# ─────────────────── 定时攒批 ───────────────────


def test_transport_drains_on_batch_interval(make_ctx, make_scalar_record):
    """After put(), the thread auto-drains and dispatches on batch_interval timeout."""
    dispatched = []
    event = threading.Event()

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            event.set()
            return True, []

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx)
    t._dispatcher = FakeDispatch()  # type: ignore
    t.start()
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1


# ─────────────────── auto_start 端到端 ───────────────────


def test_transport_integration_auto_start(make_ctx, make_scalar_record):
    """auto_start=True starts the thread and processes records normally."""
    dispatched = []
    event = threading.Event()

    class FakeDispatch:
        def __call__(self, records):
            dispatched.extend(records)
            event.set()
            return True, []

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx, auto_start=True)
    t._dispatcher = FakeDispatch()  # type: ignore
    record = make_scalar_record(step=1)
    record.num = 1
    t.put([record])
    assert event.wait(timeout=2.0)
    t.finish()
    assert len(dispatched) == 1

    # ─────────────────── consecutive failure ───────────────────


def test_transport_keeps_pending_records_and_warns_after_retry_exhaustion(make_ctx, make_scalar_record):
    """After retry exhaustion, pending records are kept and a warning is emitted."""
    attempts = []
    first_dispatch_done = threading.Event()
    success_event = threading.Event()

    class FlakyDispatch:
        def __call__(self, records):
            attempts.append([record.num for record in records])
            if len(attempts) == 1:
                first_dispatch_done.set()
                return False, records
            success_event.set()
            return True, []

    ctx = make_ctx(batch_interval=0.01)
    t = _make_transport(ctx)
    t._dispatcher = FlakyDispatch()  # type: ignore
    t.start()

    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 1
    second.num = 2

    with patch("swanlab.sdk.internal.core_python.transport.thread.UploadWarningThrottle.warn") as mock_warn:
        t.put([first])
        assert first_dispatch_done.wait(timeout=2.0)
        t.put([second])
        assert success_event.wait(timeout=2.0)
        t.finish()

    assert attempts[0] == [1]
    assert attempts[1] == [1]
    assert attempts[2] == [2]
    mock_warn.assert_called_once()
