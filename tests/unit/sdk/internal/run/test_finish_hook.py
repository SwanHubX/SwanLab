"""
@author: cunyue
@file: test_finish_hook.py
@time: 2026/3/14
@description: 测试 Run._atexit_cleanup / _excepthook / _sigint_handler 的单元行为（均 mock 依赖，不启动真实 Run）
"""

import signal
import sys
import threading
from unittest.mock import ANY, MagicMock, patch

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import ConfirmRunFinishResponse, DeliverRunFinishResponse
from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import Run, clear_run


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
    """构造一个最小化的 Run 替身"""
    mock = MagicMock(spec=Run)
    mock._state = state
    mock._init_pid = fork.current_pid()
    # alive property 依赖 fork.is_forked(_init_pid) 和 _state，手动计算
    type(mock).alive = property(lambda self: not fork.is_forked(self._init_pid) and self._state == "running")
    mock._api_lock = threading.RLock()
    return mock


def _make_real_finish_run() -> Run:
    """构造只覆盖 finish() 所需字段的 Run 实例。"""
    run = Run.__new__(Run)
    run._init_pid = fork.current_pid()
    run._state = "running"
    run._api_lock = threading.RLock()
    run.__dict__["mode"] = "online"
    run._ctx = MagicMock()
    run._components = MagicMock()
    run._probe = MagicMock()
    run._core = MagicMock()
    run._core.deliver_run_finish.return_value = DeliverRunFinishResponse(success=True, message="OK")
    run._core.confirm_run_finish.return_value = ConfirmRunFinishResponse(success=True, message="OK")
    run._callbacker = MagicMock()
    run._handle_atexit = MagicMock()
    run._sys_origin_excepthook = MagicMock()
    run._original_sigint_handler = signal.SIG_DFL
    run._is_main_thread = True
    return run


class TestAtexitCleanup:
    def test_no_op_when_not_running(self):
        """_state != 'running' 时直接返回，不调用 finish"""
        run = _make_mock_run(state="success")
        Run._handle_atexit(run)
        run.finish.assert_not_called()

    def test_calls_finish_when_running(self):
        """_state == 'running' 时应调用 finish()"""
        run = _make_mock_run(state="running")
        Run._handle_atexit(run)
        run.finish.assert_called_once()


class TestFinish:
    def test_online_finish_uses_run_with_progress(self):
        run = _make_real_finish_run()

        with (
            patch("swanlab.sdk.internal.run.greeting.goodbye"),
            patch("swanlab.sdk.internal.run.run_with_progress") as mock_run_with_progress,
            patch("swanlab.sdk.internal.run.atexit.unregister"),
            patch("swanlab.sdk.internal.run.signal.signal"),
            patch("swanlab.sdk.internal.run.console.reset"),
        ):
            mock_run_with_progress.return_value = ConfirmRunFinishResponse(success=True, message="OK")
            run.finish()

        mock_run_with_progress.assert_called_once_with(
            stats_fn=run._core.get_operation_stats,
            blocking_fn=run._core.confirm_run_finish,
        )


class TestExcepthook:
    def test_keyboard_interrupt_calls_aborted(self):
        """KeyboardInterrupt → finish(state='aborted', ...)"""
        run = _make_mock_run()
        run._sys_origin_excepthook = MagicMock()
        tp, val, tb = _make_exc_info(KeyboardInterrupt())
        Run._handle_except(run, tp, val, tb)
        run.finish.assert_called_once_with(state="aborted", error=ANY)

    def test_generic_exception_calls_crashed(self):
        """普通异常 → finish(state='crashed')，error 包含完整 traceback"""
        run = _make_mock_run()
        run._sys_origin_excepthook = MagicMock()
        tp, val, tb = _make_exc_info(RuntimeError("boom"))
        Run._handle_except(run, tp, val, tb)
        call_kwargs = run.finish.call_args.kwargs
        assert call_kwargs["state"] == "crashed"
        assert "boom" in call_kwargs["error"]

    def test_no_op_when_not_running(self):
        """_state != 'running' 时不调用 finish"""
        run = _make_mock_run(state="success")
        run._sys_origin_excepthook = MagicMock()
        tp, val, tb = _make_exc_info(RuntimeError("no run"))
        Run._handle_except(run, tp, val, tb)
        run.finish.assert_not_called()

    def test_always_calls_saved_origin_excepthook(self):
        """无论是否有活跃 Run，始终调用注册时保存的原始 hook，而非 sys.__excepthook__"""
        run = _make_mock_run(state="success")
        saved_hook = MagicMock()
        run._sys_origin_excepthook = saved_hook
        tp, val, tb = _make_exc_info(RuntimeError("test"))
        Run._handle_except(run, tp, val, tb)
        saved_hook.assert_called_once_with(tp, val, tb)

    def test_calls_saved_hook_not_builtin(self):
        """当外层框架替换了 sys.excepthook 时，调用保存的外层 hook，而非内置默认 hook"""
        run = _make_mock_run()
        outer_framework_hook = MagicMock()
        run._sys_origin_excepthook = outer_framework_hook
        tp, val, tb = _make_exc_info(RuntimeError("outer"))
        with patch("sys.__excepthook__") as mock_builtin:
            Run._handle_except(run, tp, val, tb)
        outer_framework_hook.assert_called_once_with(tp, val, tb)
        mock_builtin.assert_not_called()

    def test_internal_error_doesnt_crash(self):
        """excepthook 内部出错时不向上抛出，仍调用保存的原始 hook"""
        run = _make_mock_run()
        run.finish.side_effect = Exception("internal boom")
        saved_hook = MagicMock()
        run._sys_origin_excepthook = saved_hook
        tp, val, tb = _make_exc_info(RuntimeError("outer"))
        Run._handle_except(run, tp, val, tb)
        saved_hook.assert_called_once_with(tp, val, tb)


class TestSigintHandler:
    def test_calls_finish_aborted_when_running(self):
        """SIGINT handler 在实验运行中应调用 finish(state='aborted')，frame=None 时 error 无调用栈"""
        run = _make_mock_run()
        run._original_sigint_handler = signal.SIG_DFL
        with pytest.raises(KeyboardInterrupt):
            Run._handle_sigint(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt by user")

    def test_error_contains_stack_when_frame_provided(self):
        """frame 不为 None 时，error 应包含调用栈信息"""
        import sys

        run = _make_mock_run()
        run._original_sigint_handler = signal.SIG_IGN
        # 用当前真实 frame 模拟信号打断场景
        frame = sys._getframe()
        Run._handle_sigint(run, signal.SIGINT, frame)
        call_kwargs = run.finish.call_args.kwargs
        assert call_kwargs["state"] == "aborted"
        assert "KeyboardInterrupt by user" in call_kwargs["error"]
        # 调用栈非空时 error 应包含文件路径信息
        assert "test_finish_hook.py" in call_kwargs["error"]

    def test_no_op_when_not_running(self):
        """_state != 'running' 时不调用 finish，仍抛出 KeyboardInterrupt"""
        run = _make_mock_run(state="success")
        run._original_sigint_handler = signal.SIG_DFL
        with pytest.raises(KeyboardInterrupt):
            Run._handle_sigint(run, signal.SIGINT, None)
        run.finish.assert_not_called()

    def test_calls_original_callable_handler(self):
        """如果原始 handler 是 callable，应调用它而非 raise KeyboardInterrupt"""
        run = _make_mock_run()
        original = MagicMock()
        run._original_sigint_handler = original
        Run._handle_sigint(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt by user")
        original.assert_called_once_with(signal.SIGINT, None)

    def test_sig_ign_does_nothing(self):
        """如果原始 handler 是 SIG_IGN，应静默返回，不抛出 KeyboardInterrupt"""
        run = _make_mock_run()
        run._original_sigint_handler = signal.SIG_IGN
        # Should not raise
        Run._handle_sigint(run, signal.SIGINT, None)
        run.finish.assert_called_once_with(state="aborted", error="KeyboardInterrupt by user")


class TestSigintRegistration:
    """Run.__init__ 的 SIGINT 注册应只在主线程进行（Ray actor 等非主线程环境不注册、不崩溃）"""

    def setup_method(self):
        self._saved_excepthook = sys.excepthook
        self._original_sigint_handler = signal.getsignal(signal.SIGINT)

    def teardown_method(self):
        sys.excepthook = self._saved_excepthook
        signal.signal(signal.SIGINT, self._original_sigint_handler)
        clear_run()

    def _init_run(self) -> Run:
        """在重 mock 下执行一次真实的 Run.__init__"""
        ctx = MagicMock()
        with (
            patch("swanlab.sdk.internal.run.Components"),
            patch("swanlab.sdk.internal.run.DeliverProbeStartRequest"),
            patch("swanlab.sdk.internal.run.console"),
            patch("swanlab.sdk.internal.run.greeting"),
            patch("swanlab.sdk.internal.run.atexit"),
        ):
            return Run(ctx)

    def test_register_in_main_thread(self):
        """主线程 init：正常注册 SIGINT handler"""
        run = self._init_run()
        assert run._is_main_thread is True
        assert signal.getsignal(signal.SIGINT) == run._handle_sigint

    def test_skip_in_non_main_thread(self):
        """非主线程 init（如 Ray actor）：跳过注册且不抛 ValueError"""
        holder = {}

        def target():
            try:
                holder["run"] = self._init_run()
            except Exception as e:  # noqa: BLE001
                holder["error"] = e

        before = signal.getsignal(signal.SIGINT)
        t = threading.Thread(target=target)
        t.start()
        t.join()
        assert "error" not in holder, f"Run.__init__ raised in non-main thread: {holder.get('error')}"
        run = holder["run"]
        assert run._is_main_thread is False
        assert signal.getsignal(signal.SIGINT) == before


class TestFinishSigintRestore:
    def setup_method(self):
        self._saved_excepthook = sys.excepthook

    def teardown_method(self):
        sys.excepthook = self._saved_excepthook

    def test_restore_only_in_main_thread(self):
        """finish() 仅在主线程恢复原 handler（否则非主线程下会抛 ValueError）"""
        for is_main_thread, expected_calls in ((True, 1), (False, 0)):
            run = _make_real_finish_run()
            run._is_main_thread = is_main_thread
            with (
                patch("swanlab.sdk.internal.run.greeting.goodbye"),
                patch("swanlab.sdk.internal.run.run_with_progress") as mock_rwp,
                patch("swanlab.sdk.internal.run.atexit.unregister"),
                patch("swanlab.sdk.internal.run.signal.signal") as mock_signal,
                patch("swanlab.sdk.internal.run.console.reset"),
            ):
                mock_rwp.return_value = ConfirmRunFinishResponse(success=True, message="OK")
                run.finish()
            assert mock_signal.call_count == expected_calls
