"""RunEmitter 的背压与消费者存活保护。"""

import threading

from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus import emitter as emitter_module
from swanlab.sdk.internal.pkg import console as pkg_console


def test_emit_drops_immediately_when_consumer_is_not_running(monkeypatch):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    emitter = RunEmitter()
    emitter.set_consumer_guard(lambda: False)

    emitter.emit("event")
    emitter.emit("another")

    assert emitter.queue.empty()
    assert len(warnings) == 1
    assert "consumer is not running" in warnings[0][0]


def test_emit_enqueues_when_consumer_is_running():
    emitter = RunEmitter()
    emitter.set_consumer_guard(lambda: True)

    emitter.emit("event")

    assert emitter.queue.get_nowait() == "event"


def test_emit_blocks_while_consumer_is_running():
    emitter = RunEmitter(maxsize=1)
    emitter.set_consumer_guard(lambda: True)
    emitter.emit("first")
    produced = threading.Event()

    def produce() -> None:
        emitter.emit("second")
        produced.set()

    thread = threading.Thread(target=produce, daemon=True)
    thread.start()
    try:
        assert produced.wait(timeout=0.05) is False
        assert emitter.queue.get_nowait() == "first"
        assert produced.wait(timeout=1.0) is True
        assert emitter.queue.get_nowait() == "second"
    finally:
        thread.join(timeout=1.0)


def test_emit_unblocks_and_drops_if_consumer_exits_while_queue_full(monkeypatch):
    warnings = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: warnings.append(args))
    monkeypatch.setattr(emitter_module, "_QUEUE_PUT_POLL_INTERVAL", 0.01)
    alive = True
    emitter = RunEmitter(maxsize=1)
    emitter.set_consumer_guard(lambda: alive)
    emitter.emit("first")
    finished = threading.Event()

    def produce() -> None:
        emitter.emit("second")
        finished.set()

    thread = threading.Thread(target=produce, daemon=True)
    thread.start()
    try:
        assert finished.wait(timeout=0.05) is False
        alive = False
        assert finished.wait(timeout=1.0) is True
        assert emitter.queue.get_nowait() == "first"
        assert emitter.queue.empty()
        assert len(warnings) == 1
        assert "consumer is not running" in warnings[0][0]
    finally:
        thread.join(timeout=1.0)
