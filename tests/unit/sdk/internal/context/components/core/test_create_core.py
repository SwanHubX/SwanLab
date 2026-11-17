"""
@author: Nexisato
@file: test_create_core.py
@time: 2026/11/9
@description: create_core 的显式 opt-in 路由与二进制定位测试。
"""

from pathlib import Path
from typing import cast
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.context.components import core as core_factory
from swanlab.sdk.internal.core_go import CoreGo, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go import runtime as runtime_module
from swanlab.sdk.internal.core_go.config import DETACH_IDLE_TIMEOUT
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.settings.core import CoreSettings


def _ctx(mode: str = "online", *, attach: bool = False, idle_timeout: float = DETACH_IDLE_TIMEOUT) -> RunContext:
    ctx = MagicMock()
    ctx.config.settings.mode = mode
    ctx.config.settings.core = CoreSettings(attach=attach, idle_timeout=idle_timeout)
    return cast(RunContext, ctx)


def _use_impl(monkeypatch, impl: str) -> None:
    monkeypatch.setattr(core_factory.helper, "get_core_impl", lambda: impl)


def _write_binary(path: Path) -> Path:
    path.write_text("#!/bin/sh\n")
    return path


def test_default_uses_core_python(monkeypatch):
    _use_impl(monkeypatch, "python")

    def fail() -> Path:
        raise AssertionError("default python path must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("online"))

    assert isinstance(core, CorePython)


def test_core_opt_in_returns_core_go(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    _use_impl(monkeypatch, "go")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx("online"))

    assert isinstance(core, CoreGo)
    assert core._binary == binary  # noqa: SLF001 - 断言工厂传递的路径


@pytest.mark.parametrize("mode", ["local", "offline"])
def test_local_and_offline_use_go_core_with_binary(monkeypatch, tmp_path, mode):
    binary = _write_binary(tmp_path / "swanlab-core")
    _use_impl(monkeypatch, "go")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx(mode))

    assert isinstance(core, CoreGo)
    assert core._binary == binary  # noqa: SLF001


def test_disabled_mode_does_not_require_binary(monkeypatch):
    _use_impl(monkeypatch, "go")

    def fail() -> Path:
        raise AssertionError("disabled mode must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("disabled"))

    assert isinstance(core, CoreGo)
    assert core._binary is None  # noqa: SLF001


def test_locate_failure_propagates(monkeypatch):
    error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "missing")
    _use_impl(monkeypatch, "go")

    def fail() -> Path:
        raise error

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    with pytest.raises(CoreGoError) as caught:
        core_factory.create_core(_ctx("online"))

    assert caught.value is error


def test_create_core_passes_attach_and_idle_timeout_from_settings(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    _use_impl(monkeypatch, "go")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx("online", attach=True, idle_timeout=12.5))

    assert isinstance(core, CoreGo)
    assert core._attach is True  # noqa: SLF001
    assert core._idle_timeout == 12.5  # noqa: SLF001
