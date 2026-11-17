"""
@file: test_require.py
@description: 验证 SWANLAB_REQUIRE / swanlab.require 默认 CoreGo、no-op 告警与 core_python 回切。
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


def test_default_core_impl_is_go_and_probe_is_python():
    assert helper.get_core_impl() == "go"
    assert helper.get_probe_impl() == "python"
    assert _query_impl_in_subprocess() == "go python"


def test_env_core_python_switches_back_to_python():
    assert _query_impl_in_subprocess("core_python") == "python python"
    assert _query_impl_in_subprocess("CORE_PYTHON") == "python python"
    assert _query_impl_in_subprocess("core") == "go python"


def test_require_core_warns_noop_by_default(monkeypatch):
    warnings: List[str] = []
    monkeypatch.setattr(
        require_module.console,
        "warning",
        lambda *args, **kwargs: warnings.append(" ".join(str(a) for a in args)),
    )

    assert swanlab.require("core") == ["core"]
    assert helper.get_core_impl() == "go"
    assert len(warnings) == 1
    assert "no-op" in warnings[0]


def test_require_switches_to_core_python_and_back(monkeypatch):
    warnings: List[str] = []
    monkeypatch.setattr(
        require_module.console,
        "warning",
        lambda *args, **kwargs: warnings.append(" ".join(str(a) for a in args)),
    )

    assert swanlab.require("core_python") == ["core_python"]
    assert helper.get_core_impl() == "python"

    assert swanlab.require("CORE") == ["core"]
    assert helper.get_core_impl() == "go"
    # 从 core_python 切回 Go 是真实切换，不应提示 no-op
    assert warnings == []


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
