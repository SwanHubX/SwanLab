import threading
import time
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.timer import Timer


@pytest.fixture
def timer_manager():
    timers = []

    def _create_timer(*args, **kwargs):
        timer = Timer(*args, **kwargs)
        timers.append(timer)
        return timer

    yield _create_timer

    for timer in timers:
        try:
            timer.cancel()
            timer.join(timeout=1.0)
        except Exception:
            pass


class TestTimer:
    def test_basic_interval(self, timer_manager):
        timer = timer_manager(lambda: None, interval=0.1, immediate=False)
        timer.run()

        time.sleep(0.25)

        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count >= 2
        assert not timer.is_running

    def test_immediate_execution(self, timer_manager):
        triggered = threading.Event()

        timer = timer_manager(triggered.set, interval=10.0, immediate=True)
        timer.run()

        assert triggered.wait(0.2)

        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count == 1

    def test_dynamic_interval(self, timer_manager):
        reached_slow_interval = threading.Event()
        interval_inputs = []

        def interval_strategy(count: int) -> float:
            interval_inputs.append(count)
            if count < 2:
                return 0.1
            reached_slow_interval.set()
            return 0.5

        timer = timer_manager(lambda: None, interval=interval_strategy)
        timer.run()

        assert reached_slow_interval.wait(0.35)

        timer.cancel()
        timer.join(timeout=1.0)

        assert interval_inputs[:3] == [0, 1, 2]
        assert timer.execution_count >= 2

    def test_cancel_interrupts_sleep(self, timer_manager):
        start_time = time.monotonic()

        timer = timer_manager(lambda: None, interval=10.0)
        timer.run()

        time.sleep(0.1)

        timer.cancel()
        timer.join(timeout=1.0)

        assert time.monotonic() - start_time < 1.0
        assert not timer.is_running

    def test_error_resilience(self, timer_manager, monkeypatch):
        error_mock = MagicMock()
        monkeypatch.setattr(console, "error", error_mock)

        third_call = threading.Event()

        def failing_task() -> None:
            if timer.execution_count >= 2:
                third_call.set()
            raise ValueError("Test Error")

        timer = timer_manager(failing_task, interval=0.1)
        timer.run()

        assert third_call.wait(2.0)

        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count >= 3
        assert "Error executing task" in error_mock.call_args[0][0]

    def test_double_run_warning(self, timer_manager, monkeypatch):
        warning_mock = MagicMock()
        monkeypatch.setattr(console, "warning", warning_mock)

        entered = threading.Event()
        timer = timer_manager(entered.set, interval=1.0, immediate=True)
        timer.run()

        assert entered.wait(0.2)

        timer.run()

        warning_mock.assert_called_once_with("Timer already running")

        timer.cancel()
        timer.join(timeout=1.0)

    def test_restart_capability(self, timer_manager):
        timer = timer_manager(lambda: None, interval=0.1)

        timer.run()
        time.sleep(0.15)
        timer.cancel()
        timer.join(timeout=1.0)

        first_run_count = timer.execution_count
        assert first_run_count >= 1

        timer.run()
        time.sleep(0.15)
        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count > first_run_count

    @pytest.mark.parametrize("interval", [0, -1, False])
    def test_reject_invalid_fixed_interval(self, interval):
        with pytest.raises((TypeError, ValueError), match="Timer interval"):
            Timer(lambda: None, interval=interval)
