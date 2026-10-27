"""
@author: Nexisato
@file: test_runtime.py
@time: 2026/10/27
@description: runtime 目录创建、权限与清理测试。
"""

import os
import stat

import pytest

from swanlab.sdk.internal.core_go import process as process_module
from swanlab.sdk.internal.core_go.config import RUNTIME_DIR_PREFIX
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.process import RuntimeLease, acquire_runtime, cleanup_runtime


def test_acquire_runtime_creates_private_directory():
    lease = acquire_runtime()
    try:
        assert lease.dir.name.startswith(RUNTIME_DIR_PREFIX)
        assert lease.dir.is_dir()
        assert stat.S_IMODE(lease.dir.stat().st_mode) == 0o700
        assert lease.owner_pid == os.getpid()
        assert lease.port_file.parent == lease.dir
        assert lease.stderr_file.parent == lease.dir
    finally:
        cleanup_runtime(lease)

    assert not lease.dir.exists()


def test_cleanup_runtime_is_idempotent():
    lease = acquire_runtime()

    cleanup_runtime(lease)
    cleanup_runtime(lease)

    assert not lease.dir.exists()


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
