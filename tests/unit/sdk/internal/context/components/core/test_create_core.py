"""
@author: Nexisato
@file: test_create_core.py
@time: 2026/11/9
@description: create_core 的默认 CorePython 路由、Go core 启用与缺失回退测试。
"""

from pathlib import Path
from typing import List, Optional, cast
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.context.components import core as core_factory
from swanlab.sdk.internal.core_go import CoreGo, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go import runtime as runtime_module
from swanlab.sdk.internal.core_go.config import DETACH_IDLE_TIMEOUT
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings.core import CoreSettings
from swanlab.sdk.protocol import CoreEnum


def _ctx(
    mode: str = "online",
    *,
    detach: Optional[bool] = None,
    idle_timeout: float = DETACH_IDLE_TIMEOUT,
    parallel: str = "none",
) -> RunContext:
    ctx = MagicMock()
    ctx.config.settings.mode = mode
    ctx.config.settings.core = CoreSettings(detach=detach, idle_timeout=idle_timeout)
    ctx.config.settings.run.parallel = parallel
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


def test_default_uses_core_python_without_locating_binary(monkeypatch):
    def fail() -> Path:
        raise AssertionError("default core_python path must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("online"))

    assert helper.get_core_impl() == "python"
    assert isinstance(core, CorePython)


def test_core_go_opt_in_uses_core_go_when_binary_available(monkeypatch, tmp_path):
    helper.set_core_impl(CoreEnum.CORE)
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
    helper.set_core_impl(CoreEnum.CORE)
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = core_factory.create_core(_ctx(mode))

    assert isinstance(core, CoreGo)
    assert core._binary == binary  # noqa: SLF001


def test_disabled_mode_does_not_require_binary(monkeypatch):
    helper.set_core_impl(CoreEnum.CORE)

    def fail() -> Path:
        raise AssertionError("disabled mode must not locate the core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = core_factory.create_core(_ctx("disabled"))

    assert isinstance(core, CoreGo)
    assert core._binary is None  # noqa: SLF001


def test_missing_binary_falls_back_to_core_python_and_warns_once(monkeypatch):
    helper.set_core_impl(CoreEnum.CORE)
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
    helper.set_core_impl(CoreEnum.CORE)
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


class TestResolveDetach:
    """core.detach 未显式设置时跟随 run.parallel；显式设置优先。"""

    @staticmethod
    def _settings(detach: Optional[bool] = None, parallel: str = "none") -> Settings:
        settings = Settings(mode="online")
        settings.merge_settings({"core": {"detach": detach}, "run": {"parallel": parallel}})
        return settings

    def test_parallel_shared_shares_a_detached_core_by_default(self):
        assert core_factory.resolve_detach(self._settings(parallel="shared")) is True

    def test_parallel_none_keeps_owner_mode(self):
        assert core_factory.resolve_detach(self._settings()) is False

    def test_explicit_detach_false_overrides_parallel_shared(self):
        assert core_factory.resolve_detach(self._settings(detach=False, parallel="shared")) is False

    def test_explicit_detach_true_without_parallel(self):
        assert core_factory.resolve_detach(self._settings(detach=True)) is True

    def test_shared_parallel_reaches_core_go(self, monkeypatch, tmp_path):
        helper.set_core_impl(CoreEnum.CORE)
        binary = _write_binary(tmp_path / "swanlab-core")
        monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

        core = core_factory.create_core(_ctx("online", parallel="shared"))

        assert isinstance(core, CoreGo)
        assert core._detach is True  # noqa: SLF001 - 断言共享模式实际进入 Detached Mode


class TestIdleTimeoutDefaults:
    """DETACH_IDLE_TIMEOUT 与 CoreSettings 默认 30s 保持同步，并正确透传到 CoreGo。"""

    def test_idle_timeout_constant_matches_settings_default(self):
        assert DETACH_IDLE_TIMEOUT == 30.0
        assert CoreSettings().idle_timeout == DETACH_IDLE_TIMEOUT

    @pytest.mark.parametrize("idle_timeout", [0.0, 60.0])
    def test_parallel_shared_passes_idle_timeout_to_core_go(self, monkeypatch, tmp_path, idle_timeout):
        helper.set_core_impl(CoreEnum.CORE)
        binary = _write_binary(tmp_path / "swanlab-core")
        monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

        core = core_factory.create_core(_ctx("online", parallel="shared", idle_timeout=idle_timeout))
        assert isinstance(core, CoreGo)
        assert core._detach is True  # noqa: SLF001
        assert core._idle_timeout == idle_timeout  # noqa: SLF001

    def test_parallel_shared_uses_default_idle_timeout(self, monkeypatch, tmp_path):
        helper.set_core_impl(CoreEnum.CORE)
        binary = _write_binary(tmp_path / "swanlab-core")
        monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

        core = core_factory.create_core(_ctx("online", parallel="shared"))
        assert isinstance(core, CoreGo)
        assert core._detach is True  # noqa: SLF001
        assert core._idle_timeout == DETACH_IDLE_TIMEOUT  # noqa: SLF001


class TestDetachIgnoredForCorePython:
    """CorePython 没有独立 core 进程，detach / parallel 的共享语义必须静默忽略。"""

    @staticmethod
    def _forbid_resolve(monkeypatch) -> None:
        def fail(settings):  # noqa: ARG004 - 只用于断言未被调用
            raise AssertionError("core_python must not resolve detach")

        monkeypatch.setattr(core_factory, "resolve_detach", fail)

    def test_parallel_shared_is_silently_ignored(self, monkeypatch):
        helper.set_core_impl(CoreEnum.CORE_PYTHON)
        self._forbid_resolve(monkeypatch)
        warnings = _capture_warnings(monkeypatch)

        core = core_factory.create_core(_ctx("online", parallel="shared"))

        assert isinstance(core, CorePython)
        assert warnings == []

    def test_explicit_detach_is_silently_ignored(self, monkeypatch):
        helper.set_core_impl(CoreEnum.CORE_PYTHON)
        self._forbid_resolve(monkeypatch)
        warnings = _capture_warnings(monkeypatch)

        core = core_factory.create_core(_ctx("online", detach=True, parallel="shared"))

        assert isinstance(core, CorePython)
        assert warnings == []

    def test_missing_binary_fallback_only_warns_about_binary(self, monkeypatch):
        helper.set_core_impl(CoreEnum.CORE)
        error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "missing")
        warnings = _capture_warnings(monkeypatch)

        def fail() -> Path:
            raise error

        monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

        core = core_factory.create_core(_ctx("online", parallel="shared"))

        assert isinstance(core, CorePython)
        assert len(warnings) == 1
        assert "falling back" in warnings[0]
        assert "detach" not in warnings[0]
