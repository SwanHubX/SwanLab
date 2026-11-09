"""finish 阶段的有限重试与 Core 回收。

覆盖两层：
- ``retry_finish_call`` 只重试 Core 标注为 retryable 的失败；
- ``Run.finish`` 按 deliver → confirm → close 编排，是否在 deliver 失败后 confirm 由 Core 能力位决定。
"""

import signal
import threading
from pathlib import Path
from typing import List, Optional, cast
from unittest.mock import MagicMock, patch

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import ConfirmRunFinishResponse, DeliverRunFinishResponse
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import Run
from swanlab.sdk.internal.run.finish_retry import (
    FINISH_RETRY_ATTEMPTS,
    FINISH_RETRY_INTERVAL,
    retry_finish_call,
)
from swanlab.sdk.protocol import CoreProtocol

RUN_DIR = Path("/tmp/swanlab-finish-retry")
RUN_FILE = RUN_DIR / "run-abc.swanlab"


def _retryable(message: str = "unavailable") -> CoreGoError:
    return CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, message)


def _terminal(message: str = "rejected") -> CoreGoError:
    return CoreGoError(Stage.RPC, ReasonCode.RUN_FINISH_REJECTED, message)


class RetryCore:
    """只提供 last_error 的协议替身。"""

    def __init__(self) -> None:
        self.last_error: Optional[CoreGoError] = None


class ScriptedCore:
    """按脚本返回 deliver / confirm 响应，并记录 close 次数。"""

    confirm_after_deliver_failure = True

    def __init__(self, deliver, confirm) -> None:
        self._deliver = list(deliver)
        self._confirm = list(confirm)
        self.deliver_calls = 0
        self.confirm_calls = 0
        self.closed = 0
        self.last_error: Optional[CoreGoError] = None

    @staticmethod
    def _take(script):
        return script.pop(0) if len(script) > 1 else script[0]

    def deliver_run_finish(self, _request) -> DeliverRunFinishResponse:
        self.deliver_calls += 1
        success, message, error = self._take(self._deliver)
        self.last_error = error
        return DeliverRunFinishResponse(success=success, message=message)

    def get_operation_stats(self):
        return ConfirmRunFinishResponse(success=True, message="OK")

    def confirm_run_finish(self) -> ConfirmRunFinishResponse:
        self.confirm_calls += 1
        success, message, error = self._take(self._confirm)
        self.last_error = error
        return ConfirmRunFinishResponse(success=success, message=message)

    def close(self) -> None:
        self.closed += 1


def _make_run(core) -> Run:
    run = Run.__new__(Run)
    run._init_pid = fork.current_pid()
    run._state = "running"
    run._api_lock = threading.RLock()
    run.__dict__["mode"] = "online"
    ctx = MagicMock()
    ctx.config.settings.core.skip_store = False
    ctx.run_dir = RUN_DIR
    ctx.run_file = RUN_FILE
    run._ctx = ctx
    run._components = MagicMock()
    run._probe = MagicMock()
    run._core = core
    run._callbacker = MagicMock()
    run._handle_atexit = MagicMock()
    run._sys_origin_excepthook = MagicMock()
    run._original_sigint_handler = signal.SIG_DFL
    run._is_main_thread = True
    run._atexit_registered = True
    run._excepthook_installed = True
    run._sigint_installed = True
    run._fork_registered = True
    return run


def _finish_env(run: Run, progress: MagicMock, console: MagicMock):
    """替换 finish 依赖的外部副作用，返回 patch 上下文。"""
    return patch.multiple(
        "swanlab.sdk.internal.run",
        greeting=MagicMock(),
        run_with_progress=progress,
        atexit=MagicMock(),
        signal=MagicMock(),
        console=console,
        _current_run=run,
    )


# ── retry_finish_call ────────────────────────────────────────


def test_retry_stops_on_first_success():
    core = RetryCore()
    sleeps: List[float] = []

    result = retry_finish_call(
        cast(CoreProtocol, core),
        lambda: DeliverRunFinishResponse(success=True),
        sleep=sleeps.append,
    )

    assert result.success is True
    assert sleeps == []


def test_retry_repeats_retryable_failures_then_succeeds():
    core = RetryCore()
    calls = []
    sleeps: List[float] = []

    def call():
        calls.append(1)
        if len(calls) < 3:
            core.last_error = _retryable()
            return DeliverRunFinishResponse(success=False, message="503")
        return DeliverRunFinishResponse(success=True)

    result = retry_finish_call(cast(CoreProtocol, core), call, sleep=sleeps.append)

    assert result.success is True
    assert len(calls) == 3
    assert sleeps == [FINISH_RETRY_INTERVAL, FINISH_RETRY_INTERVAL]


def test_retry_stops_on_terminal_failure():
    core = RetryCore()
    calls = []
    sleeps: List[float] = []

    def call():
        calls.append(1)
        core.last_error = _terminal()
        return DeliverRunFinishResponse(success=False, message="rejected")

    result = retry_finish_call(cast(CoreProtocol, core), call, sleep=sleeps.append)

    assert result.success is False
    assert len(calls) == 1
    assert sleeps == []


def test_retry_is_bounded():
    core = RetryCore()
    calls = []
    sleeps: List[float] = []

    def call():
        calls.append(1)
        core.last_error = _retryable()
        return DeliverRunFinishResponse(success=False, message="503")

    result = retry_finish_call(cast(CoreProtocol, core), call, sleep=sleeps.append)

    assert result.success is False
    assert len(calls) == FINISH_RETRY_ATTEMPTS
    assert sleeps == [FINISH_RETRY_INTERVAL] * (FINISH_RETRY_ATTEMPTS - 1)


def test_retry_rejects_invalid_attempts():
    core = RetryCore()

    with pytest.raises(ValueError):
        retry_finish_call(cast(CoreProtocol, core), lambda: DeliverRunFinishResponse(success=True), attempts=0)


# ── Run.finish 编排 ──────────────────────────────────────────


def test_finish_retries_deliver_and_confirm_then_closes():
    core = ScriptedCore(
        deliver=[
            (False, "503", _retryable()),
            (True, "OK", None),
        ],
        confirm=[
            (False, "503", _retryable()),
            (False, "503", _retryable()),
            (True, "OK", None),
        ],
    )
    run = _make_run(core)
    progress = MagicMock(side_effect=lambda stats_fn, blocking_fn: blocking_fn())

    with _finish_env(run, progress, MagicMock()):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep") as mock_sleep:
            run.finish()

    assert core.deliver_calls == 2
    assert core.confirm_calls == 3
    assert core.closed == 1
    assert mock_sleep.call_count == 3


def test_finish_skips_confirm_for_core_go_when_deliver_is_terminal():
    core = ScriptedCore(
        deliver=[(False, "finish record commit failed", _terminal())],
        confirm=[(True, "OK", None)],
    )
    core.confirm_after_deliver_failure = False
    run = _make_run(core)
    progress = MagicMock()

    with _finish_env(run, progress, MagicMock()):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep"):
            run.finish()

    progress.assert_not_called()
    assert core.deliver_calls == 1
    assert core.confirm_calls == 0
    assert core.closed == 1


def test_finish_confirms_for_core_python_when_deliver_fails():
    core = ScriptedCore(
        deliver=[(False, "finish record commit failed", _terminal())],
        confirm=[(True, "OK", None)],
    )
    run = _make_run(core)
    progress = MagicMock(side_effect=lambda stats_fn, blocking_fn: blocking_fn())

    with _finish_env(run, progress, MagicMock()):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep"):
            run.finish()

    progress.assert_called_once()
    assert core.deliver_calls == 1
    assert core.confirm_calls == 1
    assert core.closed == 1


def test_finish_closes_core_when_deliver_raises():
    class RaisingCore(ScriptedCore):
        def deliver_run_finish(self, _request):
            self.deliver_calls += 1
            raise RuntimeError("boom")

    core = RaisingCore(deliver=[(True, "OK", None)], confirm=[(True, "OK", None)])
    run = _make_run(core)

    with _finish_env(run, MagicMock(), MagicMock()):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep"):
            with pytest.raises(RuntimeError, match="boom"):
                run.finish()

    assert core.closed == 1


def test_finish_reports_local_recovery_hint():
    message = "Failed to finish run, but it has been saved locally."
    core = ScriptedCore(
        deliver=[(True, "OK", None)],
        confirm=[(False, message, _retryable())],
    )
    run = _make_run(core)
    console = MagicMock()

    with _finish_env(run, MagicMock(side_effect=lambda stats_fn, blocking_fn: blocking_fn()), console):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep"):
            run.finish()

    logged = " ".join(str(arg) for arg in console.warning.call_args.args)
    assert str(RUN_FILE) in logged
    assert f"swanlab sync {RUN_DIR}" in logged


def test_finish_no_recovery_hint_with_skip_store():
    message = "Failed to finish run, but it has been saved locally."
    core = ScriptedCore(
        deliver=[(True, "OK", None)],
        confirm=[(False, message, _retryable())],
    )
    run = _make_run(core)
    run._ctx.config.settings.core.skip_store = True
    console = MagicMock()

    with _finish_env(run, MagicMock(side_effect=lambda stats_fn, blocking_fn: blocking_fn()), console):
        with patch("swanlab.sdk.internal.run.finish_retry.time.sleep"):
            run.finish()

    for call in console.warning.call_args_list:
        assert "Local copy is complete" not in " ".join(str(arg) for arg in call.args)
