"""
@author: Nexisato
@file: test_impl_sync.py
@time: 2026/11/18
@description: create_core_sync 路由、缺失二进制回退与 client.current_api_key 兜底测试。
"""

from pathlib import Path
from types import SimpleNamespace
from typing import List

from swanlab.sdk.internal import impl
from swanlab.sdk.internal.context.components import core as core_factory
from swanlab.sdk.internal.core_go import CoreGoError, CoreSyncGo, ReasonCode, Stage
from swanlab.sdk.internal.core_go import runtime as runtime_module
from swanlab.sdk.internal.core_python import client as py_client
from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.protocol import CoreEnum


def test_create_core_sync_defaults_to_core_sync_python(monkeypatch):
    def fail() -> Path:
        raise AssertionError("default python sync path must not locate core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = impl.create_core_sync()

    assert helper.get_core_impl() == "python"
    assert isinstance(core, CoreSyncPython)


def test_create_core_sync_core_go_opt_in(monkeypatch, tmp_path: Path):
    helper.set_core_impl(CoreEnum.CORE)
    binary = tmp_path / "swanlab-core"
    binary.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(runtime_module, "locate_core_binary", lambda: binary)

    core = impl.create_core_sync()

    assert helper.get_core_impl() == "go"
    assert isinstance(core, CoreSyncGo)
    assert core._binary == binary  # noqa: SLF001 - 断言工厂传递的二进制路径


def test_create_core_sync_python_opt_out_skips_locate(monkeypatch):
    helper.set_core_impl(CoreEnum.CORE_PYTHON)

    def fail() -> Path:
        raise AssertionError("python sync path must not locate core binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    core = impl.create_core_sync()

    assert isinstance(core, CoreSyncPython)


def test_create_core_sync_missing_binary_falls_back_and_warns_once(monkeypatch):
    helper.set_core_impl(CoreEnum.CORE)
    warnings: List[str] = []
    monkeypatch.setattr(
        core_factory.console,
        "warning",
        lambda *args, **kwargs: warnings.append(" ".join(str(a) for a in args)),
    )

    def fail() -> Path:
        raise CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "missing sync binary")

    monkeypatch.setattr(runtime_module, "locate_core_binary", fail)

    first = impl.create_core_sync()
    second = impl.create_core_sync()

    assert isinstance(first, CoreSyncPython)
    assert isinstance(second, CoreSyncPython)
    assert len(warnings) == 1
    assert "missing sync binary" in warnings[0]


def test_current_api_key_handles_uninitialized_and_active_client(monkeypatch):
    monkeypatch.setattr(py_client, "_default_client", None)
    assert py_client.exists() is False
    assert py_client.current_api_key() is None

    monkeypatch.setattr(py_client, "_default_client", SimpleNamespace(api_key="active-key-123"))
    assert py_client.exists() is True
    assert py_client.current_api_key() == "active-key-123"
