"""
@author: cunyue
@file: test_finish.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.cmd.finish 中各函数的单元行为（均 mock 依赖，不启动真实 Run）
"""

import sys
from unittest.mock import ANY, MagicMock

from swanlab.sdk.cmd.finish import atexit_finish, finish, swanlab_excepthook


def _make_exc_info(exc: BaseException):
    """辅助：构造 (tp, val, tb) 三元组，用于测试 excepthook"""
    try:
        raise exc
    except BaseException:
        tp, val, tb = sys.exc_info()
        assert tp is not None and val is not None
        return tp, val, tb


class TestFinishFunction:
    def test_finish_no_active_run(self, monkeypatch):
        """无活跃 Run 时，打印 error 后直接返回，不抛出异常"""
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: False)
        mock_console = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.console", mock_console)

        finish()

        mock_console.error.assert_called_once()

    def test_finish_calls_run_finish(self, monkeypatch):
        """有活跃 Run 时，应将 state 和 error 透传给 run.finish()"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", lambda: mock_run)

        finish(state="crashed", error="something went wrong")

        mock_run.finish.assert_called_once_with("crashed", "something went wrong")

    def test_finish_default_state_is_success(self, monkeypatch):
        """finish() 不传参时，默认 state 为 'success'，error 为 None"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", lambda: mock_run)

        finish()

        mock_run.finish.assert_called_once_with("success", None)


class TestAtexitFinish:
    def test_atexit_finish_no_run(self, monkeypatch):
        """无活跃 Run 时，atexit_finish 直接返回，不调用 get_run"""
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: False)
        mock_get_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", mock_get_run)

        atexit_finish()

        mock_get_run.assert_not_called()

    def test_atexit_finish_calls_run_finish(self, monkeypatch):
        """有活跃 Run 时，atexit_finish 应调用 run.finish()"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", lambda: mock_run)

        atexit_finish()

        mock_run.finish.assert_called_once()


class TestSwanlabExcepthook:
    def test_excepthook_keyboard_interrupt(self, monkeypatch):
        """KeyboardInterrupt 异常 → run.finish(state='aborted', ...)"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", lambda: mock_run)
        monkeypatch.setattr(sys, "__excepthook__", MagicMock())

        tp, val, tb = _make_exc_info(KeyboardInterrupt())
        swanlab_excepthook(tp, val, tb)

        mock_run.finish.assert_called_once_with(state="aborted", error=ANY)

    def test_excepthook_generic_exception(self, monkeypatch):
        """普通异常 → run.finish(state='crashed')，error 包含完整 traceback"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", lambda: mock_run)
        monkeypatch.setattr(sys, "__excepthook__", MagicMock())

        tp, val, tb = _make_exc_info(RuntimeError("boom"))
        swanlab_excepthook(tp, val, tb)

        call_kwargs = mock_run.finish.call_args.kwargs
        assert call_kwargs["state"] == "crashed"
        assert "boom" in call_kwargs["error"]

    def test_excepthook_no_run(self, monkeypatch):
        """无活跃 Run 时，不调用 run.finish，但不抛出异常"""
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: False)
        mock_get_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.finish.get_run", mock_get_run)
        monkeypatch.setattr(sys, "__excepthook__", MagicMock())

        tp, val, tb = _make_exc_info(RuntimeError("no run"))
        swanlab_excepthook(tp, val, tb)

        mock_get_run.assert_not_called()

    def test_excepthook_always_calls_original_hook(self, monkeypatch):
        """无论是否有活跃 Run，sys.__excepthook__ 必须被调用一次"""
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", lambda: False)
        mock_original = MagicMock()
        monkeypatch.setattr(sys, "__excepthook__", mock_original)

        tp, val, tb = _make_exc_info(RuntimeError("test"))
        swanlab_excepthook(tp, val, tb)

        mock_original.assert_called_once_with(tp, val, tb)

    def test_excepthook_internal_error_doesnt_crash(self, monkeypatch):
        """excepthook 内部逻辑出错时，不应向上抛出异常，且仍调用 sys.__excepthook__"""
        # 模拟 has_run 本身抛出异常（触发 except 块）
        monkeypatch.setattr("swanlab.sdk.cmd.finish.has_run", MagicMock(side_effect=Exception("internal boom")))
        # mock console 以隔离 console.error 的实现细节
        monkeypatch.setattr("swanlab.sdk.cmd.finish.console", MagicMock())
        mock_original = MagicMock()
        monkeypatch.setattr(sys, "__excepthook__", mock_original)

        tp, val, tb = _make_exc_info(RuntimeError("outer"))

        # 不应抛出
        swanlab_excepthook(tp, val, tb)

        # __excepthook__ 仍须被调用（finally 块保证）
        mock_original.assert_called_once_with(tp, val, tb)
