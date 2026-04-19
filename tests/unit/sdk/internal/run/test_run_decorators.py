"""
@author: cunyue
@file: test_run_decorators.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.internal.run 中的装饰器
"""

import threading

import pytest

from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import with_api


class TestWithApi:
    def test_executes_when_alive(self):
        """当 alive 为 True 且非 fork 时，方法应正常执行"""

        class MockRun:
            def __init__(self):
                self._init_pid = fork.current_pid()
                self.alive = True
                self._api_lock = threading.RLock()

            @with_api(cmd="swanlab.my_method()")
            def my_method(self, x, y):
                return x + y

        run = MockRun()
        result = run.my_method(1, 2)
        assert result == 3

    def test_raises_when_not_alive(self):
        """当 alive 为 False 时，应抛出 RuntimeError，且错误信息包含 cmd"""

        class MockRun:
            def __init__(self):
                self._init_pid = fork.current_pid()
                self.alive = False
                self._api_lock = threading.RLock()

            @with_api(cmd="swanlab.my_method()")
            def my_method(self):
                return "ok"

        run = MockRun()
        with pytest.raises(RuntimeError, match="`swanlab.my_method\\(\\)` requires an active Run"):
            run.my_method()

    def test_error_message_contains_cmd(self):

        class MockRun:
            def __init__(self):
                self._init_pid = fork.current_pid()
                self.alive = False
                self._api_lock = threading.RLock()

            @with_api(cmd="run.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="`run.log\\(\\)`"):
            run.log()

    def test_raises_when_forked(self):
        """fork 后应抛出 RuntimeError，提示用户使用 spawn"""

        class MockRun:
            def __init__(self):
                self._init_pid = 0  # 模拟 fork：PID 不匹配
                self._api_lock = threading.RLock()

            @with_api(cmd="swanlab.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="does not support fork"):
            run.log()

    def test_fork_error_takes_priority_over_not_alive(self):
        """fork 错误应优先于 not alive 错误"""

        class MockRun:
            def __init__(self):
                self._init_pid = 0
                self.alive = False
                self._api_lock = threading.RLock()

            @with_api(cmd="swanlab.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="does not support fork"):
            run.log()

    def test_must_alive_false_allows_not_alive(self):
        """must_alive=False 时，not alive 状态不抛异常"""

        class MockRun:
            def __init__(self):
                self._init_pid = fork.current_pid()
                self.alive = False
                self._api_lock = threading.RLock()

            @with_api(cmd="run.finish()", must_alive=False)
            def finish(self):
                return "ok"

        run = MockRun()
        assert run.finish() == "ok"

    def test_must_alive_false_still_raises_on_fork(self):
        """must_alive=False 时，fork 仍然抛异常"""

        class MockRun:
            def __init__(self):
                self._init_pid = 0
                self._api_lock = threading.RLock()

            @with_api(cmd="run.finish()", must_alive=False)
            def finish(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="does not support fork"):
            run.finish()

    def test_preserves_metadata(self):
        """装饰器应保留被装饰方法的 __name__ 和 __doc__"""

        class MockRun:
            @with_api(cmd="swanlab.my_method()")
            def my_method(self):
                """My docstring"""

        assert MockRun.my_method.__name__ == "my_method"
        assert MockRun.my_method.__doc__ == "My docstring"
