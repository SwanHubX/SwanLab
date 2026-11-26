"""
@file: test_throttle.py
@description: WarningThrottle 单调时钟节流器的单测。
"""

import threading

from swanlab.sdk.internal.pkg.throttle import WarningThrottle


def test_should_emit_and_reset():
    now = 100.0
    t = WarningThrottle(10.0, clock=lambda: now)

    assert t.should_emit() is True
    assert t.should_emit() is False

    now = 109.9
    assert t.should_emit() is False

    now = 110.0
    assert t.should_emit() is True
    assert t.should_emit() is False

    t.reset()
    assert t.should_emit() is True


def test_record_accumulates_and_flushes_at_window():
    now = 5.0
    t = WarningThrottle(10.0, clock=lambda: now)

    assert t.record() == 1
    assert t.record() == 0
    assert t.record(3) == 0

    now = 14.9
    assert t.record() == 0

    now = 15.0
    assert t.record() == 6
    assert t.record() == 0


def test_record_thread_safe():
    now = 1.0
    t = WarningThrottle(1000.0, clock=lambda: now)
    results = []
    lock = threading.Lock()

    def worker():
        for _ in range(50):
            val = t.record()
            with lock:
                results.append(val)

    threads = [threading.Thread(target=worker) for _ in range(4)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()

    assert sum(results) == 1
    now += 1000.0
    assert t.record() == 200
