"""
@author: cunyue
@file: test_cmd_helper.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.cmd.guard 中的装饰器
"""

import threading
import time

import pytest

from swanlab.sdk.cmd.guard import with_cmd_lock, with_run, without_run


class TestWithCmdLock:
    def test_wraps_preserves_metadata(self):
        """@wraps 应保留被装饰函数的 __name__ 和 __doc__"""

        def my_func():
            """My docstring"""

        wrapped = with_cmd_lock(my_func)

        assert wrapped.__name__ == "my_func"
        assert wrapped.__doc__ == "My docstring"

    def test_lock_serializes_concurrent_calls(self):
        """并发调用两个装饰函数时，执行顺序不得交错（严格串行化）"""
        events = []

        @with_cmd_lock
        def record_event(name):
            events.append(f"{name}:start")
            time.sleep(0.03)
            events.append(f"{name}:end")

        t1 = threading.Thread(target=record_event, args=("A",))
        t2 = threading.Thread(target=record_event, args=("B",))

        t1.start()
        time.sleep(0.005)  # 确保 t1 先抢到锁
        t2.start()

        t1.join()
        t2.join()

        # 无论谁先执行，各自的 start/end 必须紧邻，不允许交错
        assert events in [
            ["A:start", "A:end", "B:start", "B:end"],
            ["B:start", "B:end", "A:start", "A:end"],
        ]


class TestWithRun:
    def test_raises_when_no_run(self, monkeypatch):
        """无活跃 Run 时，应抛出 RuntimeError"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: False)

        @with_run("test_cmd")
        def my_func(run):
            return "ok"

        with pytest.raises(RuntimeError, match="`swanlab.test_cmd` requires an active Run"):
            my_func()

    def test_passes_run_when_active(self, monkeypatch):
        """有活跃 Run 时，装饰器不会抛出异常"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: True)

        @with_run("test_cmd")
        def my_func(x, y):
            return x + y

        result = my_func(1, 2)
        assert result == 3

    def test_preserves_metadata(self):
        """@wraps 应保留被装饰函数的 __name__ 和 __doc__"""

        @with_run("test_cmd")
        def my_func():
            """My docstring"""

        assert my_func.__name__ == "my_func"
        assert my_func.__doc__ == "My docstring"


class TestWithRunAllowFinished:
    def test_warns_when_no_run_and_allow_finished(self, monkeypatch):
        """has_run() 为 False 且 allow_finished=True 时，打印统一警告并返回 None"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: False)
        warnings = []
        monkeypatch.setattr("swanlab.sdk.cmd.guard.console.warning", lambda msg, *a, **kw: warnings.append(msg))

        @with_run("finish", allow_finished=True)
        def my_func():
            return "should not reach"

        result = my_func()
        assert result is None
        assert warnings == ["SwanLab Run has already finished or has not started."]

    def test_raises_when_no_run_and_not_allow_finished(self, monkeypatch):
        """allow_finished=False 时仍应抛异常"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: False)

        @with_run("log", allow_finished=False)
        def my_func():
            return "should not reach"

        with pytest.raises(RuntimeError, match="`swanlab.log` requires an active Run"):
            my_func()


class TestWithoutRun:
    def test_raises_when_run_exists(self, monkeypatch):
        """有活跃 Run 时，应抛出 RuntimeError"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: True)

        @without_run("test_cmd")
        def my_func():
            return "ok"

        with pytest.raises(RuntimeError, match="`swanlab.test_cmd` requires no active Run"):
            my_func()

    def test_executes_when_no_run(self, monkeypatch):
        """无活跃 Run 时，应正常执行"""
        monkeypatch.setattr("swanlab.sdk.cmd.guard.has_run", lambda: False)

        @without_run("test_cmd")
        def my_func(x, y):
            return x + y

        result = my_func(1, 2)
        assert result == 3

    def test_preserves_metadata(self):
        """@wraps 应保留被装饰函数的 __name__ 和 __doc__"""

        @without_run("test_cmd")
        def my_func():
            """My docstring"""

        assert my_func.__name__ == "my_func"
        assert my_func.__doc__ == "My docstring"
