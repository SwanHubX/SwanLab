"""
@author: Nexisato
@file: test_create_core.py
@time: 2026/11/9
@description: create_core 的默认 Go 路由、core_python 切回与缺失回退测试。
"""

from pathlib import Path
from typing import List, cast
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.context.components import core as core_factory
from swanlab.sdk.internal.core_go import CoreGo, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go import runtime as runtime_module
from swanlab.sdk.internal.core_go.config import DETACH_IDLE_TIMEOUT
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.internal.settings.core import CoreSettings
from swanlab.sdk.protocol import CoreEnum


def _ctx(mode: str = "online", *, detach: bool = False, idle_timeout: float = DETACH_IDLE_TIMEOUT) -> RunContext:
    ctx = MagicMock()
    ctx.config.settings.mode = mode
    ctx.config.settings.core = CoreSettings(detach=detach, idle_timeout=idle_timeout)
    return cast(RunContext, ctx)


def _write_binary(path: Path) -> Path:
    path.write_text("#!/bin/sh\n")
    return path


def _capture_warnings(monkeypatch) -> List[str]:
    warnings: List[str] = []
    monkeypatch.setattr(
        core_factory.console,
        "warning",
        lambda *args, **kwargs: warnings.append(" ".join(str(a) for a in args)),
    )
    return warnings


def test_default_uses_core_go_when_binary_available(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx("online"))

    assert helper.get_core_impl() == "go"
    assert isinstance(core, CoreGo)
    assert core._binary == binary  # noqa: SLF001 - 断言工厂传递的路径


def test_core_python_uses_core_python_without_locating_binary(monkeypatch):
    helper.set_core_impl(CoreEnum.CORE_PYTHON)

    def fail() -> Path:
        raise AssertionError("core_python path must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("online"))

    assert isinstance(core, CorePython)


@pytest.mark.parametrize("mode", ["local", "offline"])
def test_local_and_offline_use_go_core_with_binary(monkeypatch, tmp_path, mode):
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx(mode))

    assert isinstance(core, CoreGo)
    assert core._binary == binary  # noqa: SLF001


def test_disabled_mode_does_not_require_binary(monkeypatch):
    def fail() -> Path:
        raise AssertionError("disabled mode must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("disabled"))

    assert isinstance(core, CoreGo)
    assert core._binary is None  # noqa: SLF001


def test_missing_binary_falls_back_to_core_python_and_warns_once(monkeypatch):
    error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "missing")
    warnings = _capture_warnings(monkeypatch)

    def fail() -> Path:
        raise error

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    first = core_factory.create_core(_ctx("online"))
    second = core_factory.create_core(_ctx("offline"))

    assert isinstance(first, CorePython)
    assert isinstance(second, CorePython)
    assert len(warnings) == 1
    assert "missing" in warnings[0]


def test_create_core_passes_detach_and_idle_timeout_from_settings(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    default_core = core_factory.create_core(_ctx("online"))
    detached_core = core_factory.create_core(_ctx("online", detach=True, idle_timeout=12.5))

    assert isinstance(default_core, CoreGo)
    assert default_core._detach is False  # noqa: SLF001
    assert isinstance(detached_core, CoreGo)
    assert detached_core._detach is True  # noqa: SLF001
    assert detached_core._idle_timeout == 12.5  # noqa: SLF001


def test_create_core_uses_provided_backend(monkeypatch):
    """init 前置校验解析出的后端直接驱动创建，避免二次解析结果不一致。"""
    warnings = _capture_warnings(monkeypatch)
    reason = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "missing")
    backend = core_factory.CoreBackend("python", fallback_reason=reason)

    core = core_factory.create_core(_ctx("online"), backend=backend)

    assert isinstance(core, CorePython)
    assert len(warnings) == 1
