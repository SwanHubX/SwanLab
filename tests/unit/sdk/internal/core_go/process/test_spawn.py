"""
@author: Nexisato
@file: test_spawn.py
@time: 2026/10/27
@description: Popen 参数、环境变量、早死分类与 port-file 等待测试。
"""

import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go import process as process_module
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.process import (
    acquire_runtime,
    cleanup_runtime,
    reap_process,
    spawn_core,
    stderr_tail,
    wait_endpoint,
)


@pytest.fixture(autouse=True)
def isolated_runtime_root(tmp_path, monkeypatch):
    """隔离测试用 runtime 根目录，避免写入真实 per-user runtime 目录。"""
    if sys.platform == "win32":
        monkeypatch.setattr(process_module, "_runtime_root", lambda: tmp_path)
        yield
        return

    root = Path(tempfile.mkdtemp(prefix="swl-t-", dir="/tmp"))
    monkeypatch.setattr(process_module, "_runtime_root", lambda: root)
    try:
        yield
    finally:
        shutil.rmtree(root, ignore_errors=True)


@pytest.fixture(autouse=True)
def fast_reap(monkeypatch):
    """这些用例只验证等待逻辑，回收超时压缩到亚秒级。"""
    monkeypatch.setattr(process_module, "WAIT_TIMEOUT", 0.2)
    monkeypatch.setattr(process_module, "TERMINATE_TIMEOUT", 0.2)
    monkeypatch.setattr(process_module, "KILL_TIMEOUT", 2.0)


def _wait_for(predicate, timeout: float = 2.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    return False


def _write_script(path: Path, body: str) -> Path:
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(0o755)
    return path


@pytest.mark.skipif(sys.platform == "win32", reason="fake core is a POSIX shell script")
def test_spawn_core_passes_argv_and_sanitizes_env(tmp_path, monkeypatch):
    argv_file = tmp_path / "argv.txt"
    env_file = tmp_path / "env.txt"
    script = _write_script(
        tmp_path / "fake-core",
        f'printf "%s\\n" "$@" > "{argv_file}"\nenv > "{env_file}"\necho "stderr-line" >&2\nsleep 30\n',
    )
    monkeypatch.setenv("SWANLAB_CORE_LISTEN", "unix:/tmp/other.sock")
    monkeypatch.setenv("SWANLAB_CORE_PARENT_PID", "1")

    lease = acquire_runtime()
    proc = None
    try:
        proc = spawn_core(script, lease)

        assert _wait_for(argv_file.exists), "fake core did not start"
        assert _wait_for(lambda: "stderr-line" in stderr_tail(lease))
        argv = argv_file.read_text().splitlines()
        assert argv == [
            "--port-filename",
            str(lease.port_file),
            "--parent-pid",
            str(os.getpid()),
        ]
        env = env_file.read_text()
        assert "SWANLAB_CORE_LISTEN=" not in env
        assert "SWANLAB_CORE_PARENT_PID=" not in env
    finally:
        reap_process(proc)
        cleanup_runtime(lease)


def test_spawn_core_resolves_relative_binary(tmp_path, monkeypatch):
    """回归：子进程 cwd 是 runtime 目录，相对路径必须提前解析为绝对路径。"""
    binary = tmp_path / "fake-core"
    binary.write_text("")
    binary.chmod(0o755)

    captured = {}

    class FakePopen:
        def __init__(self, *args, **kwargs):
            captured.update(kwargs)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(process_module.subprocess, "Popen", FakePopen)

    lease = acquire_runtime()
    try:
        spawn_core(Path("fake-core"), lease)
    finally:
        cleanup_runtime(lease)

    assert Path(captured["args"][0]) == binary.resolve()


def test_spawn_core_rejects_missing_binary(tmp_path):
    lease = acquire_runtime()
    try:
        with pytest.raises(CoreGoError) as caught:
            spawn_core(tmp_path / "missing-core", lease)

        assert caught.value.stage is Stage.SPAWN
        assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    finally:
        cleanup_runtime(lease)


def test_spawn_core_rejects_non_executable_binary(tmp_path):
    binary = tmp_path / "fake-core"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o644)
    lease = acquire_runtime()
    try:
        with pytest.raises(CoreGoError) as caught:
            spawn_core(binary, lease)

        assert caught.value.reason is ReasonCode.SPAWN_FAILED
    finally:
        cleanup_runtime(lease)


@pytest.mark.parametrize(
    "exit_code, expected",
    [
        (3, ReasonCode.ENDPOINT_BUSY),
        (2, ReasonCode.USAGE_ERROR),
        (1, ReasonCode.EARLY_EXIT),
    ],
)
def test_wait_endpoint_classifies_early_exit(exit_code, expected):
    lease = acquire_runtime()
    proc = subprocess.Popen([sys.executable, "-c", f"import sys; sys.exit({exit_code})"])
    try:
        with pytest.raises(CoreGoError) as caught:
            wait_endpoint(proc, lease, timeout=2.0)

        assert caught.value.stage is Stage.PORT_FILE
        assert caught.value.reason is expected
    finally:
        reap_process(proc)
        cleanup_runtime(lease)


def test_wait_endpoint_times_out_without_port_file():
    lease = acquire_runtime()
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        with pytest.raises(CoreGoError) as caught:
            wait_endpoint(proc, lease, timeout=0.2)

        assert caught.value.reason is ReasonCode.PORT_FILE_TIMEOUT
    finally:
        reap_process(proc)
        cleanup_runtime(lease)


def test_wait_endpoint_rejects_pid_mismatch_immediately():
    lease = acquire_runtime()
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        # Windows 的 Go core 发布 sock=<port>，unix 路径只在 POSIX 下合法。
        endpoint = "sock=12345" if sys.platform == "win32" else f"unix={lease.dir / 'core.sock'}"
        lease.port_file.write_bytes(f"protocol=1\npid={proc.pid + 1}\n{endpoint}\nEOF\n".encode())

        started = time.monotonic()
        with pytest.raises(CoreGoError) as caught:
            wait_endpoint(proc, lease, timeout=5.0)

        assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH
        assert time.monotonic() - started < 1.0
    finally:
        reap_process(proc)
        cleanup_runtime(lease)


def test_wait_endpoint_ignores_partial_port_file():
    lease = acquire_runtime()
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        lease.port_file.write_bytes(b"protocol=1\npid=")

        with pytest.raises(CoreGoError) as caught:
            wait_endpoint(proc, lease, timeout=0.2)

        assert caught.value.reason is ReasonCode.PORT_FILE_TIMEOUT
    finally:
        reap_process(proc)
        cleanup_runtime(lease)


@pytest.mark.parametrize(
    "idle_timeout, expected_tail",
    [
        (2.0, ["--detach", "--idle-timeout", "2s"]),
        (0.25, ["--detach", "--idle-timeout", "250ms"]),
        (0.0, ["--detach"]),
    ],
)
def test_spawn_core_detach_passes_flags_without_parent_pid(tmp_path, monkeypatch, idle_timeout, expected_tail):
    binary = tmp_path / "fake-core"
    binary.write_text("")
    binary.chmod(0o755)

    captured = {}

    class FakePopen:
        def __init__(self, *args, **kwargs):
            captured.update(kwargs)

    monkeypatch.setattr(process_module.subprocess, "Popen", FakePopen)

    lease = acquire_runtime()
    try:
        spawn_core(binary, lease, detach=True, idle_timeout=idle_timeout)
    finally:
        cleanup_runtime(lease)

    argv = captured["args"]
    assert "--parent-pid" not in argv
    assert argv[1:3] == ["--port-filename", str(lease.port_file)]
    assert argv[3:] == expected_tail


def test_spawn_core_rejects_negative_idle_timeout(tmp_path):
    binary = tmp_path / "fake-core"
    binary.write_text("")
    binary.chmod(0o755)
    lease = acquire_runtime()
    try:
        with pytest.raises(ValueError):
            spawn_core(binary, lease, detach=True, idle_timeout=-1.0)
    finally:
        cleanup_runtime(lease)
