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
        counter = 0
        lock = threading.Lock()

        def task() -> None:
            nonlocal counter
            with lock:
                counter += 1

        timer = timer_manager(task, interval=0.1, immediate=False)
        timer.run()

        time.sleep(0.25)

        timer.cancel()
        timer.join(timeout=1.0)

        assert counter >= 2
        assert timer.execution_count >= 2
        assert not timer.is_running

    def test_immediate_execution(self, timer_manager):
        triggered = threading.Event()
        counter = 0

        def task() -> None:
            nonlocal counter
            counter += 1
            triggered.set()

        timer = timer_manager(task, interval=10.0, immediate=True)
        timer.run()

        assert triggered.wait(0.2)

        time.sleep(0.05)

        timer.cancel()
        timer.join(timeout=1.0)

        assert counter == 1
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

        assert 0 in interval_inputs
        assert 1 in interval_inputs
        assert 2 in interval_inputs
        assert timer.execution_count >= 2

    def test_cancel_interrupts_sleep(self, timer_manager):
        start_time = time.monotonic()

        timer = timer_manager(lambda: None, interval=10.0)
        timer.run()

        time.sleep(0.1)

        timer.cancel()
        timer.join(timeout=1.0)

        duration = time.monotonic() - start_time

        assert duration < 1.0
        assert not timer.is_running

    def test_error_resilience(self, timer_manager, monkeypatch):
        error_mock = MagicMock()
        monkeypatch.setattr(console, "error", error_mock)

        mock_task = MagicMock(side_effect=ValueError("Test Error"))

        timer = timer_manager(mock_task, interval=0.1)
        timer.run()

        time.sleep(0.35)

        timer.cancel()
        timer.join(timeout=1.0)

        assert mock_task.call_count >= 3
        error_mock.assert_called()
        assert "Error executing task" in error_mock.call_args[0][0]
        assert timer.execution_count >= 3

    def test_double_run_warning(self, timer_manager, monkeypatch):
        warning_mock = MagicMock()
        entered = threading.Event()
        released = threading.Event()

        monkeypatch.setattr(console, "warning", warning_mock)

        def task() -> None:
            entered.set()
            released.wait(0.5)

        timer = timer_manager(task, interval=1.0, immediate=True)
        timer.run()

        assert entered.wait(0.2)

        timer.run()

        warning_mock.assert_called_once_with("Timer already running")

        timer.cancel()
        released.set()
        timer.join(timeout=1.0)

    def test_restart_capability(self, timer_manager):
        counter = 0

        def task() -> None:
            nonlocal counter
            counter += 1

        timer = timer_manager(task, interval=0.1)

        timer.run()
        time.sleep(0.15)
        timer.cancel()
        timer.join(timeout=1.0)

        first_run_count = counter
        assert first_run_count >= 1

        timer.run()
        time.sleep(0.15)
        timer.cancel()
        timer.join(timeout=1.0)

        assert counter > first_run_count

    @pytest.mark.parametrize("interval", [0, -1, False])
    def test_reject_invalid_fixed_interval(self, interval):
        with pytest.raises((TypeError, ValueError), match="Timer interval"):
            Timer(lambda: None, interval=interval)
