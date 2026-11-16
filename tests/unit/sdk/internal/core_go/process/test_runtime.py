"""
@author: Nexisato
@file: test_runtime.py
@time: 2026/10/27
@description: per-user runtime 目录创建、权限校验与清理测试。
"""

import os
import shutil
import stat
import sys
import tempfile
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go import process as process_module
from swanlab.sdk.internal.core_go.config import RUNTIME_DIR_MODE, RUNTIME_DIR_PREFIX
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.process import (
    RuntimeLease,
    acquire_runtime,
    cleanup_runtime,
    default_runtime_dir,
)


@pytest.fixture
def short_root(tmp_path):
    if sys.platform == "win32":
        yield tmp_path
        return
    root = Path(tempfile.mkdtemp(prefix="swl-rt-", dir="/tmp"))
    try:
        yield root
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_acquire_runtime_creates_and_reuses_per_user_directory(short_root, monkeypatch):
    monkeypatch.setattr(process_module, "_runtime_root", lambda: short_root)

    first = acquire_runtime()
    second = acquire_runtime()
    try:
        assert first.dir == default_runtime_dir(short_root)
        assert first.dir == second.dir
        assert first.dir.name.startswith(RUNTIME_DIR_PREFIX)
        assert first.dir.is_dir()
        if sys.platform != "win32":
            assert stat.S_IMODE(first.dir.stat().st_mode) == RUNTIME_DIR_MODE
        assert first.owner_pid == os.getpid()
        assert first.port_file.parent == first.dir
        assert first.stderr_file.parent == first.dir
    finally:
        cleanup_runtime(first)

    # 共享 per-user runtime 目录不得被单个客户端清理删除。
    assert first.dir.exists()


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX permissions and symlinks")
def test_acquire_runtime_rejects_bad_perm_and_symlink(short_root, monkeypatch):
    monkeypatch.setattr(process_module, "_runtime_root", lambda: short_root)
    target = default_runtime_dir(short_root)

    other = short_root / "other"
    other.mkdir()
    target.symlink_to(other)
    with pytest.raises(CoreGoError) as symlink_err:
        acquire_runtime()
    assert symlink_err.value.stage is Stage.RUNTIME
    assert symlink_err.value.reason is ReasonCode.RUNTIME_UNAVAILABLE

    target.unlink()
    target.mkdir(mode=0o755)
    os.chmod(target, 0o755)
    with pytest.raises(CoreGoError) as perm_err:
        acquire_runtime()
    assert perm_err.value.reason is ReasonCode.RUNTIME_UNAVAILABLE
    # 对已存在的不合规目录失败时不得擅自删除。
    assert target.exists()


def test_cleanup_runtime_cleans_legacy_mkdtemp_directory(tmp_path):
    legacy_dir = tmp_path / f"{RUNTIME_DIR_PREFIX}legacy-rand"
    legacy_dir.mkdir()
    lease = RuntimeLease(
        dir=legacy_dir,
        port_file=legacy_dir / "core.port",
        stderr_file=legacy_dir / "core.stderr.log",
        owner_pid=os.getpid(),
    )

    cleanup_runtime(lease)
    cleanup_runtime(lease)

    assert not legacy_dir.exists()


def test_cleanup_runtime_skips_foreign_owner(tmp_path):
    runtime_dir = tmp_path / f"{RUNTIME_DIR_PREFIX}foreign"
    runtime_dir.mkdir()
    lease = RuntimeLease(
        dir=runtime_dir,
        port_file=runtime_dir / "core.port",
        stderr_file=runtime_dir / "core.stderr.log",
        owner_pid=os.getpid() + 1,
    )

    cleanup_runtime(lease)

    assert runtime_dir.exists()


def test_cleanup_runtime_skips_unprefixed_directory(tmp_path):
    runtime_dir = tmp_path / "not-a-runtime"
    runtime_dir.mkdir()
    lease = RuntimeLease(
        dir=runtime_dir,
        port_file=runtime_dir / "core.port",
        stderr_file=runtime_dir / "core.stderr.log",
        owner_pid=os.getpid(),
    )

    cleanup_runtime(lease)

    assert runtime_dir.exists()


@pytest.mark.skipif(sys.platform == "win32", reason="UDS path limit is POSIX-only")
def test_acquire_runtime_rejects_oversized_socket_path(tmp_path, monkeypatch):
    root = tmp_path / "runtime-root"
    root.mkdir()
    monkeypatch.setattr(process_module, "_runtime_root", lambda: root)
    monkeypatch.setattr(process_module, "MAX_UNIX_PATH_BYTES", 8)

    with pytest.raises(CoreGoError) as caught:
        acquire_runtime()

    assert caught.value.stage is Stage.RUNTIME
    assert caught.value.reason is ReasonCode.RUNTIME_UNAVAILABLE
    assert list(root.iterdir()) == []


@pytest.mark.parametrize(
    "env, expected",
    [
        ({"USERNAME": r"CORP\Alice.Smith"}, "CORP_Alice.Smith"),
        ({"USERNAME": "张三"}, "default"),
        ({"USERNAME": "alice_张三"}, "alice"),
        ({"USERNAME": "...___"}, "default"),
        ({"USERNAME": "  ", "USER": "bob", "LOGNAME": "carol"}, "bob"),
        ({}, "default"),
    ],
)
def test_windows_user_segment_matches_go_contract(env, expected):
    assert process_module._windows_user_segment(env) == expected
