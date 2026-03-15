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
                self._state = "running"

            @with_run
            def my_method(self, x, y):
                return x + y

        run = MockRun()
        result = run.my_method(1, 2)
        assert result == 3

    def test_raises_when_state_is_not_running(self):
        """当 _state 不为 'running' 时，应抛出 RuntimeError"""

        class MockRun:
            def __init__(self):
                self._state = "finished"

            @with_run
            def my_method(self):
                return "ok"

        run = MockRun()
        with pytest.raises(RuntimeError, match="`swanlab.run` requires an active SwanLabRun"):
            run.my_method()

    def test_preserves_metadata(self):
        """装饰器应保留被装饰方法的 __name__ 和 __doc__"""

        class MockRun:
            @with_run
            def my_method(self):
                """My docstring"""

        assert MockRun.my_method.__name__ == "my_method"
        assert MockRun.my_method.__doc__ == "My docstring"
