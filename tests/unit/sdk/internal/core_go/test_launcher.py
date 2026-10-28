"""
@author: Nexisato
@file: test_launcher.py
@time: 2026/10/27
@description: CoreLauncher 的启动编排、失败回收、幂等关闭与 fork guard。
"""

import os
from pathlib import Path
from typing import Any, Optional, cast

import pytest

from swanlab.sdk.internal.core_go import launcher as launcher_module
from swanlab.sdk.internal.core_go.discovery import Endpoint
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.process import RuntimeLease


class FakeProcess:
    pid = 4242


class FakeChannel:
    last: Optional["FakeChannel"] = None

    def __init__(self, endpoint: Endpoint) -> None:
        self.endpoint = endpoint
        self.closed = 0

    @classmethod
    def open(cls, endpoint: Endpoint, timeout: float = 0.0) -> "FakeChannel":
        cls.last = cls(endpoint)
        return cls.last

    def close(self) -> None:
        self.closed += 1


@pytest.fixture
def launcher_env(tmp_path, monkeypatch):
    binary = tmp_path / "swanlab-core"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    lease = RuntimeLease(
        dir=tmp_path / "runtime",
        port_file=tmp_path / "runtime" / "core.port",
        stderr_file=tmp_path / "runtime" / "core.stderr.log",
        owner_pid=os.getpid(),
    )
    lease.dir.mkdir()
    endpoint = Endpoint(kind="unix", pid=FakeProcess.pid, unix_path=str(lease.dir / "core.sock"))
    calls = {"acquire": 0, "spawn": 0, "wait": 0, "reap": 0, "cleanup": 0}

    def fake_acquire() -> RuntimeLease:
        calls["acquire"] += 1
        return lease

    def fake_spawn(binary_path: Path, runtime: RuntimeLease) -> FakeProcess:
        calls["spawn"] += 1
        return FakeProcess()

    def fake_wait(proc, runtime: RuntimeLease, timeout: float = 0.0) -> Endpoint:
        calls["wait"] += 1
        return endpoint

    def fake_reap(proc) -> int:
        calls["reap"] += 1
        return 0

    def fake_cleanup(runtime: RuntimeLease) -> None:
        calls["cleanup"] += 1

    monkeypatch.setattr(launcher_module, "acquire_runtime", fake_acquire)
    monkeypatch.setattr(launcher_module, "spawn_core", fake_spawn)
    monkeypatch.setattr(launcher_module, "wait_endpoint", fake_wait)
    monkeypatch.setattr(launcher_module, "reap_process", fake_reap)
    monkeypatch.setattr(launcher_module, "cleanup_runtime", fake_cleanup)
    monkeypatch.setattr(launcher_module, "Channel", FakeChannel)

    yield binary, lease, endpoint, calls

    FakeChannel.last = None


def test_start_acquires_runtime_and_opens_channel(launcher_env):
    binary, lease, endpoint, calls = launcher_env
    launcher = CoreLauncher(binary)

    assert launcher.start() == endpoint
    assert calls == {"acquire": 1, "spawn": 1, "wait": 1, "reap": 0, "cleanup": 0}
    assert launcher.endpoint == endpoint
    assert launcher.channel is FakeChannel.last


def test_start_is_idempotent(launcher_env):
    binary, _, endpoint, calls = launcher_env
    launcher = CoreLauncher(binary)

    assert launcher.start() == endpoint
    assert launcher.start() == endpoint
    assert calls["acquire"] == 1
    assert calls["spawn"] == 1


def test_injected_lease_skips_acquire(launcher_env):
    binary, lease, _, calls = launcher_env
    launcher = CoreLauncher(binary, lease=lease)

    launcher.start()

    assert calls["acquire"] == 0
    assert calls["spawn"] == 1


def test_close_cleans_channel_process_and_runtime_once(launcher_env):
    binary, _, _, calls = launcher_env
    launcher = CoreLauncher(binary)
    launcher.start()
    channel = cast(Any, launcher.channel)

    launcher.close()
    launcher.close()

    assert channel.closed == 1
    assert calls["reap"] == 1
    assert calls["cleanup"] == 1


def test_abort_is_idempotent(launcher_env):
    binary, _, _, calls = launcher_env
    launcher = CoreLauncher(binary)
    launcher.start()
    channel = cast(Any, launcher.channel)

    launcher.abort()
    launcher.abort()

    assert channel.closed == 1
    assert calls["reap"] == 1
    assert calls["cleanup"] == 1


def test_start_missing_binary_does_not_acquire_runtime(launcher_env, tmp_path):
    _, _, _, calls = launcher_env
    launcher = CoreLauncher(tmp_path / "missing-core")

    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.stage is Stage.SPAWN
    assert caught.value.reason is ReasonCode.BINARY_NOT_FOUND
    assert calls["acquire"] == 0
    assert calls["cleanup"] == 0


def test_start_failure_after_lease_cleans_up(launcher_env, monkeypatch):
    binary, _, _, calls = launcher_env

    def fail_spawn(binary_path: Path, runtime: RuntimeLease):
        raise CoreGoError(Stage.SPAWN, ReasonCode.SPAWN_FAILED, "spawn boom")

    monkeypatch.setattr(launcher_module, "spawn_core", fail_spawn)
    launcher = CoreLauncher(binary)

    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.reason is ReasonCode.SPAWN_FAILED
    assert calls["reap"] == 0
    assert calls["cleanup"] == 1
    with pytest.raises(CoreGoError) as closed:
        _ = launcher.endpoint
    assert closed.value.reason is ReasonCode.CLOSED


def test_start_channel_failure_reaps_process_and_runtime(launcher_env, monkeypatch):
    binary, _, _, calls = launcher_env

    def fail_open(endpoint: Endpoint, timeout: float = 0.0):
        raise CoreGoError(Stage.CHANNEL, ReasonCode.CHANNEL_TIMEOUT, "channel boom")

    monkeypatch.setattr(FakeChannel, "open", fail_open)
    launcher = CoreLauncher(binary)

    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.reason is ReasonCode.CHANNEL_TIMEOUT
    assert calls["reap"] == 1
    assert calls["cleanup"] == 1


def test_fork_child_cannot_start_or_reap_parent_core(launcher_env):
    binary, _, _, calls = launcher_env
    launcher = CoreLauncher(binary)
    launcher._owner_pid = os.getpid() + 1

    with pytest.raises(CoreGoError) as caught:
        launcher.start()
    assert caught.value.reason is ReasonCode.FORK_UNSUPPORTED

    launcher.close()
    launcher.abort()

    assert calls["spawn"] == 0
    assert calls["reap"] == 0
    assert calls["cleanup"] == 0


def test_fork_child_does_not_close_running_parent_core(launcher_env):
    binary, _, _, calls = launcher_env
    launcher = CoreLauncher(binary)
    launcher.start()
    channel = cast(Any, launcher.channel)

    launcher._owner_pid = os.getpid() + 1
    launcher.close()
    launcher.abort()

    assert channel.closed == 0
    assert calls["reap"] == 0
    assert calls["cleanup"] == 0


def test_properties_report_not_started_and_closed(launcher_env):
    binary, _, _, _ = launcher_env
    launcher = CoreLauncher(binary)

    with pytest.raises(CoreGoError) as not_started:
        _ = launcher.endpoint
    assert not_started.value.reason is ReasonCode.NOT_STARTED

    launcher.start()
    launcher.close()

    with pytest.raises(CoreGoError) as closed:
        _ = launcher.channel
    assert closed.value.reason is ReasonCode.CLOSED
