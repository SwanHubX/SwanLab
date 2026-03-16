"""
@author: cunyue
@file: test_finish_hook.py
@time: 2026/3/14
@description: 测试 SwanLabRun._atexit_cleanup / _excepthook / _sigint_handler 的单元行为（均 mock 依赖，不启动真实 Run）
"""

import signal
import sys
import threading
from unittest.mock import ANY, MagicMock, patch

import pytest

from swanlab.sdk.internal.run import SwanLabRun


def _make_exc_info(exc: BaseException):
    """辅助：构造 (tp, val, tb) 三元组"""
    # noinspection PyBroadException
    try:
        raise exc
    except BaseException:
        tp, val, tb = sys.exc_info()
        assert tp is not None and val is not None
        return tp, val, tb


def _make_mock_run(state: str = "running") -> MagicMock:
    """构造一个最小化的 SwanLabRun 替身"""
    mock = MagicMock(spec=SwanLabRun)
    mock._state = state
    mock._api_lock = threading.RLock()
    return mock


class TestAtexitCleanup:
    def test_no_op_when_not_running(self):
        """_state != 'running' 时直接返回，不调用 finish"""
        run = _make_mock_run(state="success")
        SwanLabRun._atexit_cleanup(run)
        run.finish.assert_not_called()

    def test_calls_finish_when_running(self):
        """_state == 'running' 时应调用 finish()"""
        run = _make_mock_run(state="running")
        SwanLabRun._atexit_cleanup(run)
        run.finish.assert_called_once()


class TestExcepthook:
    def test_keyboard_interrupt_calls_aborted(self):
        """KeyboardInterrupt → finish(state='aborted', ...)"""
        run = _make_mock_run()
        with patch("sys.__excepthook__"):
            tp, val, tb = _make_exc_info(KeyboardInterrupt())
            SwanLabRun._excepthook(run, tp, val, tb)
        run.finish.assert_called_once_with(state="aborted", error=ANY)

    def test_generic_exception_calls_crashed(self):
        """普通异常 → finish(state='crashed')，error 包含完整 traceback"""
        run = _make_mock_run()
        with patch("sys.__excepthook__"):
            tp, val, tb = _make_exc_info(RuntimeError("boom"))
            SwanLabRun._excepthook(run, tp, val, tb)
        call_kwargs = run.finish.call_args.kwargs
        assert call_kwargs["state"] == "crashed"
        assert "boom" in call_kwargs["error"]

    def test_no_op_when_not_running(self):
        """_state != 'running' 时不调用 finish"""
        run = _make_mock_run(state="success")
        with patch("sys.__excepthook__"):
            tp, val, tb = _make_exc_info(RuntimeError("no run"))
            SwanLabRun._excepthook(run, tp, val, tb)
        run.finish.assert_not_called()

    def test_always_calls_original_excepthook(self):
        """无论是否有活跃 Run，sys.__excepthook__ 必须被调用一次"""
        run = _make_mock_run(state="success")
        with patch("sys.__excepthook__") as mock_original:
            tp, val, tb = _make_exc_info(RuntimeError("test"))
            SwanLabRun._excepthook(run, tp, val, tb)
            mock_original.assert_called_once_with(tp, val, tb)

    def test_internal_error_doesnt_crash(self):
        """excepthook 内部出错时不向上抛出，仍调用 sys.__excepthook__"""
        run = _make_mock_run()
        run.finish.side_effect = Exception("internal boom")
        with patch("sys.__excepthook__") as mock_original:
            tp, val, tb = _make_exc_info(RuntimeError("outer"))
            SwanLabRun._excepthook(run, tp, val, tb)
            mock_original.assert_called_once_with(tp, val, tb)


class TestSigintHandler:
    def test_calls_finish_aborted_when_running(self):
        """SIGINT handler 在实验运行中应调用 finish(state='aborted')"""
        run = _make_mock_run()
        run._original_sigint_handler = signal.SIG_DFL
        with pytest.raises(KeyboardInterrupt):
            SwanLabRun._sigint_handler(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt")

    def test_no_op_when_not_running(self):
        """_state != 'running' 时不调用 finish，仍抛出 KeyboardInterrupt"""
        run = _make_mock_run(state="success")
        run._original_sigint_handler = signal.SIG_DFL
        with pytest.raises(KeyboardInterrupt):
            SwanLabRun._sigint_handler(run, signal.SIGINT, None)
        run.finish.assert_not_called()

    def test_calls_original_callable_handler(self):
        """如果原始 handler 是 callable，应调用它而非 raise KeyboardInterrupt"""
        run = _make_mock_run()
        original = MagicMock()
        run._original_sigint_handler = original
        SwanLabRun._sigint_handler(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt")
        original.assert_called_once_with(signal.SIGINT, None)

    def test_sig_ign_does_nothing(self):
        """如果原始 handler 是 SIG_IGN，应静默返回，不抛出 KeyboardInterrupt"""
        run = _make_mock_run()
        run._original_sigint_handler = signal.SIG_IGN
        # Should not raise
        SwanLabRun._sigint_handler(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt")
