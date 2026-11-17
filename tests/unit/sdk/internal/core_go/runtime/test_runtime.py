"""
@author: Nexisato
@file: test_runtime.py
@time: 2026/10/27
@description: per-user runtime 目录创建、权限校验与清理测试。
"""

import os
import stat
import sys
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go import runtime as runtime_module
from swanlab.sdk.internal.core_go.config import (
    CORE_BINARY_ENV,
    RUNTIME_DIR_MODE,
    RUNTIME_DIR_PREFIX,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.runtime import (
    RuntimeLease,
    acquire_runtime,
    cleanup_runtime,
    default_runtime_dir,
    locate_core_binary,
)


def test_acquire_runtime_creates_and_reuses_per_user_directory(short_tmp_path, monkeypatch):
    monkeypatch.setattr(runtime_module, "_runtime_root", lambda: short_tmp_path)

    first = acquire_runtime()
    second = acquire_runtime()
    try:
        assert first.dir == default_runtime_dir(short_tmp_path)
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


def test_acquire_owner_runtime_isolates_directories_and_cleans_up(short_tmp_path, monkeypatch):
    monkeypatch.setattr(runtime_module, "_runtime_root", lambda: short_tmp_path)

    shared = acquire_runtime(shared=True)
    owner_a = acquire_runtime(shared=False)
    owner_b = acquire_runtime(shared=False)
    try:
        assert owner_a.dir != shared.dir
        assert owner_b.dir != shared.dir
        assert owner_a.dir != owner_b.dir
        assert owner_a.dir.name.startswith(f"{RUNTIME_DIR_PREFIX}own-")
        assert owner_a.dir.is_dir()
        if sys.platform != "win32":
            assert stat.S_IMODE(owner_a.dir.stat().st_mode) == RUNTIME_DIR_MODE
    finally:
        cleanup_runtime(owner_a)
        cleanup_runtime(owner_b)
        cleanup_runtime(shared)

    assert not owner_a.dir.exists()
    assert not owner_b.dir.exists()
    assert shared.dir.exists()


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX permissions and symlinks")
def test_acquire_runtime_rejects_bad_perm_and_symlink(short_tmp_path, monkeypatch):
    monkeypatch.setattr(runtime_module, "_runtime_root", lambda: short_tmp_path)
    target = default_runtime_dir(short_tmp_path)

    other = short_tmp_path / "other"
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
    monkeypatch.setattr(runtime_module, "_runtime_root", lambda: root)
    monkeypatch.setattr(runtime_module, "MAX_UNIX_PATH_BYTES", 8)

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
    assert runtime_module._windows_user_segment(env) == expected


def _write_binary(path: Path) -> Path:
    path.write_text("#!/bin/sh\n")
    return path


def test_locate_env_override_is_preferred(monkeypatch, tmp_path):
    binary = _write_binary(tmp_path / "custom-core")
    monkeypatch.setenv(CORE_BINARY_ENV, str(binary))

    assert locate_core_binary() == binary.resolve()


def test_locate_env_override_missing_raises(monkeypatch, tmp_path):
    missing = tmp_path / "missing-core"
    monkeypatch.setenv(CORE_BINARY_ENV, str(missing))

    with pytest.raises(CoreGoError) as caught:
        locate_core_binary()

    assert caught.value.stage is Stage.SPAWN
    assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    assert CORE_BINARY_ENV in caught.value.message


def test_locate_blank_env_override_falls_back_to_packaged(monkeypatch, tmp_path):
    packaged = _write_binary(tmp_path / "swanlab-core")
    monkeypatch.setenv(CORE_BINARY_ENV, "   ")
    monkeypatch.setattr(runtime_module, "_PACKAGED_BINARY", packaged)

    assert locate_core_binary() == packaged.resolve()


def test_locate_missing_packaged_binary_raises(monkeypatch, tmp_path):
    monkeypatch.delenv(CORE_BINARY_ENV, raising=False)
    monkeypatch.setattr(runtime_module, "_PACKAGED_BINARY", tmp_path / "swanlab-core")

    with pytest.raises(CoreGoError) as caught:
        locate_core_binary()

    assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    assert "platform wheel" in caught.value.message


def test_packaged_binary_points_into_swanlab_bin():
    packaged = runtime_module._PACKAGED_BINARY

    assert packaged.parent.name == "bin"
    assert packaged.parent.parent.name == "swanlab"
    assert packaged.name.startswith("swanlab-core")
