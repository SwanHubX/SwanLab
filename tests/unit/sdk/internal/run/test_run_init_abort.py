"""Run.__init__ 失败回滚与 fork 子进程 Core 脱离。

关键契约：
- init 任一步失败都必须停组件/probe、abort Core、撤销 hook 与全局 Run；
- 失败不得留下运行中的 Run，也不得覆盖原始异常；
- fork 子进程只拿到 _DetachedCoreGo，不操作父 Core。
"""

import signal
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.core_go.core import CoreGo, _DetachedCoreGo
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.run import Run, has_run
from swanlab.sdk.protocol import CoreProtocol


def _make_ctx() -> MagicMock:
    ctx = MagicMock()
    settings = ctx.config.settings
    settings.mode = "online"
    settings.project.name = "project"
    settings.run.id = "run-id"
    settings.core.skip_store = False
    settings.model_copy.return_value = MagicMock()
    settings.to_probe_proto.return_value = MagicMock()
    ctx.run_dir = Path("/tmp/swanlab-run-init-abort")
    ctx.global_system_step = 0
    return ctx


def _init_patches(ctx: MagicMock, components: MagicMock, **overrides):
    """装配 Run.__init__ 的外部依赖，返回 patch 上下文。"""
    patcher = patch.multiple(
        "swanlab.sdk.internal.run",
        Components=MagicMock(return_value=components),
        DeliverProbeStartRequest=MagicMock(),
        console=MagicMock(),
        greeting=overrides.pop("greeting", MagicMock()),
        atexit=MagicMock(),
        **overrides,
    )
    return patcher


def test_run_init_failure_aborts_core_and_cleans_hooks():
    ctx = _make_ctx()
    components = MagicMock()
    greeting = MagicMock()
    greeting.welcome.side_effect = RuntimeError("welcome failed")
    original_excepthook = sys.excepthook

    with (
        _init_patches(ctx, components, greeting=greeting),
        patch("swanlab.sdk.internal.run.fork.register"),
        patch("swanlab.sdk.internal.run.fork.unregister") as mock_unregister,
        patch("swanlab.sdk.internal.run.signal.signal") as mock_signal,
        pytest.raises(RuntimeError, match="welcome failed"),
    ):
        Run(ctx)

    ctx.core.abort.assert_called_once()
    ctx.probe.deliver_probe_finish.assert_called_once()
    assert components.start.called
    assert components.stop.called
    mock_unregister.assert_called_once()
    mock_signal.assert_called_with(signal.SIGINT, signal.getsignal(signal.SIGINT))
    assert sys.excepthook is original_excepthook
    assert not has_run()


def test_run_init_failure_before_components_start_does_not_stop_components():
    ctx = _make_ctx()
    components = MagicMock()
    ctx.callbacker.on_run_initialized.side_effect = RuntimeError("callback failed")

    with (
        _init_patches(ctx, components),
        patch("swanlab.sdk.internal.run.fork.register"),
        patch("swanlab.sdk.internal.run.fork.unregister"),
        patch("swanlab.sdk.internal.run.signal.signal"),
        pytest.raises(RuntimeError, match="callback failed"),
    ):
        Run(ctx)

    ctx.core.abort.assert_called_once()
    components.start.assert_not_called()
    components.stop.assert_not_called()
    ctx.probe.deliver_probe_finish.assert_not_called()
    assert not has_run()


def test_run_init_failure_during_components_start_still_stops_components():
    ctx = _make_ctx()
    components = MagicMock()
    components.start.side_effect = RuntimeError("start failed")

    with (
        _init_patches(ctx, components),
        patch("swanlab.sdk.internal.run.fork.register"),
        patch("swanlab.sdk.internal.run.fork.unregister"),
        patch("swanlab.sdk.internal.run.signal.signal"),
        pytest.raises(RuntimeError, match="start failed"),
    ):
        Run(ctx)

    ctx.core.abort.assert_called_once()
    components.stop.assert_called_once()
    ctx.probe.deliver_probe_finish.assert_not_called()


def test_run_init_failure_during_probe_start_finishes_probe():
    ctx = _make_ctx()
    components = MagicMock()
    ctx.probe.deliver_probe_start.side_effect = RuntimeError("probe failed")

    with (
        _init_patches(ctx, components),
        patch("swanlab.sdk.internal.run.fork.register"),
        patch("swanlab.sdk.internal.run.fork.unregister"),
        patch("swanlab.sdk.internal.run.signal.signal"),
        pytest.raises(RuntimeError, match="probe failed"),
    ):
        Run(ctx)

    ctx.core.abort.assert_called_once()
    components.stop.assert_called_once()
    ctx.probe.deliver_probe_finish.assert_called_once()


def test_abort_init_is_idempotent():
    run = Run.__new__(Run)
    run._state = "running"
    run._core = MagicMock(spec=CoreProtocol)
    run._probe = MagicMock()
    run._components = MagicMock()
    run._atexit_registered = False
    run._excepthook_installed = False
    run._sigint_installed = False
    run._fork_registered = False
    run._probe_started = True
    run._components_started = True

    with patch("swanlab.sdk.internal.run.console.reset"):
        run._abort_init(RuntimeError("boom"))
        run._abort_init(RuntimeError("boom"))

    run._core.abort.assert_called_once()  # type: ignore[attr-defined]
    run._probe.deliver_probe_finish.assert_called_once()
    run._components.stop.assert_called_once()


def test_handle_fork_replaces_core_go_with_detached():
    run = Run.__new__(Run)
    parent = CoreGo("online")
    run._core = parent

    Run._handle_fork(run)

    assert isinstance(run._core, _DetachedCoreGo)
    assert run._core is not parent
    assert parent.last_error is None


def test_handle_fork_keeps_core_python():
    run = Run.__new__(Run)
    core = CorePython("online")
    run._core = core

    Run._handle_fork(run)

    assert run._core is core


def test_handle_fork_swallows_detach_errors():
    run = Run.__new__(Run)
    core = MagicMock()
    core.fork.side_effect = RuntimeError("detach failed")
    run._core = core

    Run._handle_fork(run)

    assert run._core is core
