"""
@author: Nexisato
@file: test_locate.py
@time: 2026/11/9
@description: swanlab-core 二进制定位优先级与错误信息测试。
"""

from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go import locate as locate_module
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.locate import CORE_BINARY_ENV, locate_core_binary


def _write_binary(path: Path) -> Path:
    path.write_text("#!/bin/sh\n")
    return path


def test_env_override_is_preferred(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "custom-core")
    monkeypatch.setenv(CORE_BINARY_ENV, str(binary))

    assert locate_core_binary() == binary.resolve()


def test_env_override_missing_raises(monkeypatch, tmp_path):
    missing = tmp_path / "missing-core"
    monkeypatch.setenv(CORE_BINARY_ENV, str(missing))

    with pytest.raises(CoreGoError) as caught:
        locate_core_binary()

    assert caught.value.stage is Stage.SPAWN
    assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    assert CORE_BINARY_ENV in caught.value.message


def test_empty_env_override_falls_back_to_packaged(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setenv(CORE_BINARY_ENV, "   ")
    monkeypatch.setattr(locate_module, "_packaged_binary", lambda: binary)

    assert locate_core_binary() == binary


def test_packaged_binary_is_fallback(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.delenv(CORE_BINARY_ENV, raising=False)
    monkeypatch.setattr(locate_module, "_packaged_binary", lambda: binary)

    assert locate_core_binary() == binary


def test_missing_binary_raises_with_hint(monkeypatch, tmp_path):
    missing = tmp_path / "swanlab-core"
    monkeypatch.delenv(CORE_BINARY_ENV, raising=False)
    monkeypatch.setattr(locate_module, "_packaged_binary", lambda: missing)

    with pytest.raises(CoreGoError) as caught:
        locate_core_binary()

    assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    assert "bundled binary" in caught.value.message
    assert CORE_BINARY_ENV in caught.value.message


def test_packaged_path_is_inside_swanlab_bin(monkeypatch):
    monkeypatch.setattr(locate_module, "_IS_WINDOWS", False)

    packaged = locate_module._packaged_binary()

    assert packaged.parent.name == "bin"
    assert packaged.parent.parent.name == "swanlab"
    assert packaged.name == "swanlab-core"


def test_windows_binary_uses_exe_suffix(monkeypatch):
    monkeypatch.setattr(locate_module, "_IS_WINDOWS", True)

    assert locate_module._binary_name() == "swanlab-core.exe"
