"""
@author: cunyue
@file: test_decorators.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.internal.run 中的装饰器
"""

import pytest

from swanlab.sdk.internal.run import with_run


class TestWithRun:
    def test_executes_when_state_is_running(self):
        """当 _state 为 'running' 时，方法应正常执行"""

        class MockRun:
            def __init__(self):
                self._forked = False
                self.alive = True

            @with_run(cmd="swanlab.my_method()")
            def my_method(self, x, y):
                return x + y

        run = MockRun()
        result = run.my_method(1, 2)
        assert result == 3

    def test_raises_when_state_is_not_running(self):
        """当 _state 不为 'running' 时，应抛出 RuntimeError，且错误信息包含 cmd"""

        class MockRun:
            def __init__(self):
                self._forked = False
                self.alive = False

            @with_run(cmd="swanlab.my_method()")
            def my_method(self):
                return "ok"

        run = MockRun()
        with pytest.raises(RuntimeError, match="`swanlab.my_method\\(\\)` requires an active Run"):
            run.my_method()

    def test_error_message_contains_cmd(self):
        """错误信息中应包含传入的 cmd 字符串"""

        class MockRun:
            def __init__(self):
                self._forked = False
                self.alive = False

            @with_run(cmd="run.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="`run.log\\(\\)`"):
            run.log()

    def test_raises_when_forked(self):
        """fork 后应抛出 RuntimeError，提示用户使用 spawn"""

        class MockRun:
            def __init__(self):
                self._forked = True
                self.alive = False

            @with_run(cmd="swanlab.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="does not support fork"):
            run.log()

    def test_fork_error_takes_priority_over_not_running(self):
        """fork 错误应优先于 not running 错误"""

        class MockRun:
            def __init__(self):
                self._forked = True
                self.alive = False

            @with_run(cmd="swanlab.log()")
            def log(self):
                pass

        run = MockRun()
        with pytest.raises(RuntimeError, match="does not support fork"):
            run.log()

    def test_preserves_metadata(self):
        """装饰器应保留被装饰方法的 __name__ 和 __doc__"""

        class MockRun:
            @with_run(cmd="swanlab.my_method()")
            def my_method(self):
                """My docstring"""

        assert MockRun.my_method.__name__ == "my_method"
        assert MockRun.my_method.__doc__ == "My docstring"
