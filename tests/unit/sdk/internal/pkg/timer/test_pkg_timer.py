import threading
import time
from typing import cast
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
        """验证固定间隔下任务被周期执行，且 execution_count 正确递增。"""
        task_triggered = threading.Event()

        def task():
            if timer.execution_count >= 1:
                task_triggered.set()

        timer = timer_manager(task, interval=0.1, immediate=False)
        timer.run()

        assert task_triggered.wait(timeout=2.0), "Task was not executed twice"

        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count >= 2
        assert not timer.is_running

    def test_immediate_execution(self, timer_manager):
        """验证 immediate=True 时任务在启动后立即执行一次，不计入后续等待周期。"""
        task_triggered = threading.Event()

        timer = timer_manager(task_triggered.set, interval=10.0, immediate=True)
        timer.run()

        assert task_triggered.wait(2.0)

        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count == 1

    def test_dynamic_interval(self, timer_manager):
        """验证动态间隔策略接收已执行次数并据此调整间隔，策略入参递增。"""
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

        assert reached_slow_interval.wait(2.0)

        timer.cancel()
        timer.join(timeout=1.0)

        assert interval_inputs[:3] == [0, 1, 2]
        assert timer.execution_count >= 2

    def test_cancel_interrupts_sleep(self, timer_manager):
        """验证 cancel() 能立即唤醒正在等待间隔的线程，不会等到间隔结束。"""
        start_time = time.monotonic()

        timer = timer_manager(lambda: None, interval=10.0)
        timer.run()

        time.sleep(0.1)

        timer.cancel()
        timer.join(timeout=1.0)

        assert time.monotonic() - start_time < 1.0
        assert not timer.is_running

    def test_error_resilience(self, timer_manager, monkeypatch):
        """验证任务抛异常不影响后续轮次执行，异常信息通过 console.trace 输出。"""
        trace_mock = MagicMock()
        monkeypatch.setattr(console, "trace", trace_mock)

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
        assert "Error executing task" in trace_mock.call_args[0][0]

    def test_double_run_warning(self, timer_manager, monkeypatch):
        """验证运行中再次调用 run() 只告警不启动重复线程。"""
        warning_mock = MagicMock()
        monkeypatch.setattr(console, "warning", warning_mock)

        task_triggered = threading.Event()
        timer = timer_manager(task_triggered.set, interval=1.0, immediate=True)
        timer.run()

        assert task_triggered.wait(2.0)

        timer.run()

        warning_mock.assert_called_once_with("Timer already running")

        timer.cancel()
        timer.join(timeout=1.0)

    def test_restart_capability(self, timer_manager):
        """验证 cancel + join 后可以重新 run()，execution_count 持续累加。"""
        task_triggered = threading.Event()

        timer = timer_manager(task_triggered.set, interval=0.1)

        timer.run()
        assert task_triggered.wait(timeout=2.0), "First run: task was never executed"
        timer.cancel()
        timer.join(timeout=1.0)

        first_run_count = timer.execution_count
        assert first_run_count >= 1

        task_triggered.clear()
        timer.run()
        assert task_triggered.wait(timeout=2.0), "Second run: task was never executed"
        timer.cancel()
        timer.join(timeout=1.0)

        assert timer.execution_count > first_run_count

    @pytest.mark.parametrize("interval", [0, -1, False])
    def test_reject_invalid_fixed_interval(self, interval):
        """验证构造时传入 0、负数或 bool 作为间隔会抛出 TypeError/ValueError。"""
        with pytest.raises(cast(tuple[type[Exception], ...], (TypeError, ValueError)), match="Timer interval"):
            Timer(lambda: None, interval=interval)

    def test_interval_strategy_runtime_error(self, timer_manager, monkeypatch):
        """验证动态间隔策略在运行时抛异常时，循环终止并输出错误日志。"""
        trace_mock = MagicMock()
        monkeypatch.setattr(console, "trace", trace_mock)

        call_count = 0

        def failing_strategy(count: int) -> float:
            nonlocal call_count
            call_count += 1
            if call_count >= 2:
                raise RuntimeError("strategy exploded")
            return 0.05

        timer = timer_manager(lambda: None, interval=failing_strategy)
        timer.run()

        # 等待循环因策略异常而退出
        deadline = time.monotonic() + 3.0
        while timer.is_running and time.monotonic() < deadline:
            time.sleep(0.05)

        assert not timer.is_running, "Timer should have stopped after interval strategy error"
        assert "Timer interval strategy error" in trace_mock.call_args[0][0]

    def test_cancel_before_start_has_no_effect(self, timer_manager):
        """验证在启动前调用 cancel() 不会阻止后续启动（start() 会重置 stop_event）。"""
        task_triggered = threading.Event()

        timer = timer_manager(task_triggered.set, interval=1.0, immediate=True)
        timer.cancel()
        timer.run()

        assert task_triggered.wait(timeout=2.0), "Timer should start normally despite prior cancel"
        timer.cancel()
        timer.join(timeout=1.0)
        assert timer.execution_count >= 1
