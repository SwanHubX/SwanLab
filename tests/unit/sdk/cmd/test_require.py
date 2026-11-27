"""
@file: test_require.py
@description: 验证 SWANLAB_REQUIRE / swanlab.require 默认 CorePython、静默幂等切换与 core 启用。
"""

import os
import subprocess
import sys
from typing import List, Optional

import pytest

import swanlab
from swanlab.sdk.cmd import require as require_module
from swanlab.sdk.internal.pkg import helper


def _query_impl_in_subprocess(env_value: Optional[str] = None) -> str:
    env = {k: v for k, v in os.environ.items() if k != "SWANLAB_REQUIRE"}
    if env_value is not None:
        env["SWANLAB_REQUIRE"] = env_value
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "from swanlab.sdk.internal.pkg import helper; print(helper.get_core_impl(), helper.get_probe_impl())",
        ],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return proc.stdout.strip()


@pytest.fixture
def captured_warnings(monkeypatch) -> List[str]:
    captured: List[str] = []
    monkeypatch.setattr(
        require_module.console,
        "warning",
        lambda *args, **kwargs: captured.append(" ".join(str(a) for a in args)),
    )
    return captured


def test_default_core_impl_is_python_and_probe_is_python():
    assert helper.get_core_impl() == "python"
    assert helper.get_probe_impl() == "python"
    assert _query_impl_in_subprocess() == "python python"


def test_env_core_switches_to_go_and_core_python_stays_python():
    assert _query_impl_in_subprocess("core") == "go python"
    assert _query_impl_in_subprocess("CORE") == "go python"
    assert _query_impl_in_subprocess("core_python") == "python python"
    assert _query_impl_in_subprocess("CORE_PYTHON") == "python python"


def test_require_is_quietly_idempotent(captured_warnings):
    assert swanlab.require("core_python") == ["core_python"]
    assert helper.get_core_impl() == "python"

    assert swanlab.require("core") == ["core"]
    assert helper.get_core_impl() == "go"

    assert swanlab.require("core") == ["core"]
    assert helper.get_core_impl() == "go"
    assert captured_warnings == []


def test_require_switches_to_core_and_back(captured_warnings):
    assert helper.get_core_impl() == "python"

    assert swanlab.require("CORE") == ["core"]
    assert helper.get_core_impl() == "go"

    assert swanlab.require("core_python") == ["core_python"]
    assert helper.get_core_impl() == "python"
    assert captured_warnings == []


def test_require_unknown_token_raises():
    with pytest.raises(ValueError, match="Unknown requirement"):
        swanlab.require("nonexistent")


def test_require_raises_when_run_active():
    run = swanlab.init(mode="disabled")
    try:
        with pytest.raises(RuntimeError, match="requires no active Run"):
            swanlab.require("core_python")
    finally:
        run.finish()
