"""
@author: Nexisato
@file: test_launcher.py
@time: 2026/10/27
@description: CoreLauncher 的启动编排、失败回收、幂等关闭、二进制定位与 fork guard。
"""

import os
import sys
import threading
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Optional, cast

import pytest

from swanlab.sdk.internal.core_go import launcher as launcher_module
from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.runtime import Endpoint, RuntimeLease, acquire_start_lock


class FakeProcess:
    pid = 4242

    def __init__(self) -> None:
        self.returncode: Optional[int] = None

    def poll(self) -> Optional[int]:
        return self.returncode


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
    lease.dir.mkdir(mode=0o700)
    if sys.platform != "win32":
        os.chmod(lease.dir, 0o700)
    endpoint = Endpoint(kind="unix", pid=FakeProcess.pid, unix_path=str(lease.dir / "core.sock"))
    calls: dict[str, Any] = {
        "acquire": 0,
        "acquire_shared": [],
        "lock": 0,
        "discover": 0,
        "spawn": 0,
        "wait": 0,
        "reap": 0,
        "cleanup": 0,
    }

    def fake_acquire(*, shared: bool = True) -> RuntimeLease:
        calls["acquire"] += 1
        calls["acquire_shared"].append(shared)
        return lease

    def fake_lock(runtime_dir: Path):
        calls["lock"] += 1
        return nullcontext()

    def fake_discover(port_file: Path, runtime_dir: Path) -> Optional[Endpoint]:
        calls["discover"] += 1
        return None

    def fake_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs) -> FakeProcess:
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
    monkeypatch.setattr(launcher_module, "acquire_start_lock", fake_lock)
    monkeypatch.setattr(launcher_module, "discover_endpoint", fake_discover)
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
    assert calls == {
        "acquire": 1,
        "acquire_shared": [False],
        "lock": 0,
        "discover": 0,
        "spawn": 1,
        "wait": 1,
        "reap": 0,
        "cleanup": 0,
    }
    assert launcher.owns_process is True
    assert launcher.endpoint == endpoint
    assert launcher.channel is FakeChannel.last


def test_process_exited_reports_owned_child_exit(launcher_env):
    binary, _, _, _ = launcher_env
    launcher = CoreLauncher(binary)

    assert launcher.process_exited is False
    launcher.start()
    assert launcher.process_exited is False

    cast(FakeProcess, launcher._proc).returncode = 1  # noqa: SLF001
    assert launcher.process_exited is True


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


@pytest.mark.parametrize("action", ["close", "abort"])
def test_close_and_abort_clean_channel_process_and_runtime_once(launcher_env, action):
    binary, _, _, calls = launcher_env
    launcher = CoreLauncher(binary)
    launcher.start()
    channel = cast(Any, launcher.channel)

    getattr(launcher, action)()
    getattr(launcher, action)()

    assert channel.closed == 1
    assert calls["reap"] == 1
    assert calls["cleanup"] == 1
    assert launcher.exit_code == 0


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

    def fail_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs):
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


def test_detach_reuses_existing_endpoint_without_lock_or_spawn(launcher_env, monkeypatch):
    binary, _, endpoint, calls = launcher_env
    monkeypatch.setattr(launcher_module, "discover_endpoint", lambda *_: endpoint)

    launcher = CoreLauncher(binary, detach=True)
    assert launcher.start() == endpoint
    assert calls["acquire_shared"] == [True]
    assert calls["lock"] == 0
    assert calls["spawn"] == 0
    assert launcher.owns_process is False

    channel = cast(Any, launcher.channel)
    launcher.close()
    assert channel.closed == 1
    assert calls["reap"] == 0
    assert calls["cleanup"] == 0
    assert launcher.exit_code is None


def test_detach_double_check_under_lock_reuses_winner_endpoint(launcher_env, monkeypatch):
    binary, _, endpoint, calls = launcher_env
    answers = iter([None, endpoint])
    monkeypatch.setattr(launcher_module, "discover_endpoint", lambda *_: next(answers))

    launcher = CoreLauncher(binary, detach=True)
    assert launcher.start() == endpoint
    assert calls["lock"] == 1
    assert calls["spawn"] == 0
    assert launcher.owns_process is False

    launcher.close()
    assert calls["reap"] == 0
    assert calls["cleanup"] == 0


def test_detach_converges_on_endpoint_busy(launcher_env, monkeypatch):
    binary, _, endpoint, calls = launcher_env
    answers = iter([None, None, endpoint])
    monkeypatch.setattr(launcher_module, "discover_endpoint", lambda *_: next(answers))

    def busy_wait(proc, runtime: RuntimeLease, timeout: float = 0.0):
        raise CoreGoError(Stage.PORT_FILE, ReasonCode.ENDPOINT_BUSY, "busy")

    monkeypatch.setattr(launcher_module, "wait_endpoint", busy_wait)

    launcher = CoreLauncher(binary, detach=True)
    assert launcher.start() == endpoint
    assert calls["spawn"] == 1
    assert calls["reap"] == 1
    assert launcher.owns_process is False

    launcher.close()
    assert calls["reap"] == 1
    assert calls["cleanup"] == 0


def test_detach_reuse_channel_failure_falls_back_to_spawn(launcher_env, monkeypatch):
    """复用刚好撞上旧 Core 空闲关闭导致 Channel 失败时，应回退到持锁启动新实例。"""
    binary, _, endpoint, calls = launcher_env
    answers = iter([endpoint, None])
    monkeypatch.setattr(launcher_module, "discover_endpoint", lambda *_: next(answers))

    open_calls = 0

    def flaky_open(ep: Endpoint, timeout: float = 0.0):
        nonlocal open_calls
        open_calls += 1
        if open_calls == 1:
            raise CoreGoError(Stage.CHANNEL, ReasonCode.GRPC_UNAVAILABLE, "dying core")
        return FakeChannel(ep)

    monkeypatch.setattr(launcher_module.Channel, "open", flaky_open)

    launcher = CoreLauncher(binary, detach=True)
    assert launcher.start() == endpoint
    assert open_calls == 2
    assert calls["lock"] == 1
    assert calls["spawn"] == 1
    assert launcher.owns_process is True

    launcher.close()


def test_detach_reuse_missing_grpcio_fails_without_spawn(launcher_env, monkeypatch):
    """缺少 grpcio 属于本机环境问题：不删 port-file、不拉新进程，直接抛出。"""
    binary, lease, endpoint, calls = launcher_env
    lease.port_file.write_text("4242 unix:/tmp/core.sock\n")
    monkeypatch.setattr(launcher_module, "discover_endpoint", lambda *_: endpoint)

    def missing_grpc(ep: Endpoint, timeout: float = 0.0):
        raise CoreGoError(Stage.CHANNEL, ReasonCode.GRPC_LOCAL_ERROR, "grpcio is required")

    monkeypatch.setattr(launcher_module.Channel, "open", missing_grpc)

    launcher = CoreLauncher(binary, detach=True)
    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.reason is ReasonCode.GRPC_LOCAL_ERROR
    assert calls["lock"] == 0
    assert calls["spawn"] == 0
    assert lease.port_file.exists()


def test_detach_raises_endpoint_busy_when_fallback_discover_misses(launcher_env, monkeypatch):
    binary, _, _, calls = launcher_env

    def busy_wait(proc, runtime: RuntimeLease, timeout: float = 0.0):
        raise CoreGoError(Stage.PORT_FILE, ReasonCode.ENDPOINT_BUSY, "busy")

    monkeypatch.setattr(launcher_module, "wait_endpoint", busy_wait)

    launcher = CoreLauncher(binary, detach=True)
    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.reason is ReasonCode.ENDPOINT_BUSY
    assert calls["reap"] == 1
    assert calls["cleanup"] == 0


def test_concurrent_detach_launchers_with_real_start_lock_spawn_once(launcher_env, monkeypatch):
    binary, lease, endpoint, calls = launcher_env
    published: dict[str, Optional[Endpoint]] = {"endpoint": None}

    def stateful_discover(port_file: Path, runtime_dir: Path) -> Optional[Endpoint]:
        calls["discover"] += 1
        return published["endpoint"]

    def publish_wait(proc, runtime: RuntimeLease, timeout: float = 0.0) -> Endpoint:
        calls["wait"] += 1
        published["endpoint"] = endpoint
        return endpoint

    monkeypatch.setattr(launcher_module, "acquire_start_lock", acquire_start_lock)
    monkeypatch.setattr(launcher_module, "discover_endpoint", stateful_discover)
    monkeypatch.setattr(launcher_module, "wait_endpoint", publish_wait)

    launchers = [CoreLauncher(binary, detach=True) for _ in range(4)]
    barrier = threading.Barrier(len(launchers))
    results = []
    errors = []
    guard = threading.Lock()

    def worker(item: CoreLauncher):
        try:
            barrier.wait()
            ep = item.start()
            with guard:
                results.append(ep)
        except BaseException as exc:  # pragma: no cover
            with guard:
                errors.append(exc)

    threads = [threading.Thread(target=worker, args=(item,)) for item in launchers]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert errors == []
    assert results == [endpoint] * 4
    assert calls["spawn"] == 1
    assert sum(1 for item in launchers if item.owns_process) == 1

    for item in launchers:
        item.close()
    assert calls["reap"] == 0


def test_detach_clears_stale_port_file_before_spawning(launcher_env, monkeypatch):
    binary, lease, endpoint, calls = launcher_env
    lease.port_file.write_bytes(b"protocol=1\npid=99999\nunix=/tmp/old.sock\nEOF\n")
    saw_port_file_at_spawn = {}

    def checking_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs):
        saw_port_file_at_spawn["exists"] = runtime.port_file.exists()
        calls["spawn"] += 1
        return FakeProcess()

    monkeypatch.setattr(launcher_module, "spawn_core", checking_spawn)

    launcher = CoreLauncher(binary, detach=True)
    assert launcher.start() == endpoint
    assert saw_port_file_at_spawn == {"exists": False}
    launcher.close()


@pytest.mark.parametrize("first_mode", ["attach", "owner"])
def test_owner_runs_are_isolated_from_each_other_and_shared_attach(short_tmp_path, monkeypatch, first_mode):
    from swanlab.sdk.internal.core_go import runtime as runtime_module

    root = short_tmp_path
    try:
        binary = root / "swanlab-core"
        binary.write_text("#!/bin/sh\n")
        binary.chmod(0o755)
        monkeypatch.setattr(runtime_module, "_runtime_root", lambda: root)
        monkeypatch.setattr(launcher_module, "acquire_runtime", runtime_module.acquire_runtime)
        monkeypatch.setattr(launcher_module, "acquire_start_lock", runtime_module.acquire_start_lock)
        monkeypatch.setattr(launcher_module, "cleanup_runtime", runtime_module.cleanup_runtime)
        monkeypatch.setattr(launcher_module, "Channel", FakeChannel)

        spawned: list[tuple[Path, dict[str, Any]]] = []
        reaped: list[int] = []
        published: dict[Path, Endpoint] = {}

        class TaggedProc:
            def __init__(self, pid: int, runtime_dir: Path) -> None:
                self.pid = pid
                self.runtime_dir = runtime_dir

        def fake_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs) -> TaggedProc:
            proc = TaggedProc(pid=5000 + len(spawned), runtime_dir=runtime.dir)
            spawned.append((runtime.dir, kwargs))
            return proc

        def fake_wait(proc: TaggedProc, runtime: RuntimeLease, timeout: float = 0.0) -> Endpoint:
            ep = Endpoint(kind="tcp", pid=proc.pid, tcp_port=30000 + proc.pid)
            runtime.port_file.write_bytes(f"protocol=1\npid={proc.pid}\nsock={ep.tcp_port}\nEOF\n".encode("ascii"))
            published[runtime.dir] = ep
            return ep

        def fake_discover(port_file: Path, runtime_dir: Path) -> Optional[Endpoint]:
            return published.get(runtime_dir)

        def fake_reap(proc: TaggedProc) -> int:
            reaped.append(proc.pid)
            return 0

        monkeypatch.setattr(launcher_module, "spawn_core", fake_spawn)
        monkeypatch.setattr(launcher_module, "wait_endpoint", fake_wait)
        monkeypatch.setattr(launcher_module, "discover_endpoint", fake_discover)
        monkeypatch.setattr(launcher_module, "reap_process", fake_reap)

        attach = CoreLauncher(binary, detach=True, idle_timeout=2.0)
        attach_peer = CoreLauncher(binary, detach=True, idle_timeout=2.0)
        owner_a = CoreLauncher(binary, detach=False)
        owner_b = CoreLauncher(binary, detach=False)

        if first_mode == "attach":
            attach_ep = attach.start()
            attach_peer_ep = attach_peer.start()
            owner_a_ep = owner_a.start()
            owner_b_ep = owner_b.start()
        else:
            owner_a_ep = owner_a.start()
            owner_b_ep = owner_b.start()
            attach_ep = attach.start()
            attach_peer_ep = attach_peer.start()

        # Attach 客户端共享 1 个 Core；每个 Owner 客户端各自独占 1 个 Core 和隔离目录。
        assert len(spawned) == 3
        assert attach_peer_ep == attach_ep
        assert owner_a_ep != owner_b_ep
        assert owner_a_ep != attach_ep
        assert owner_b_ep != attach_ep
        assert owner_a._lease is not None and owner_b._lease is not None and attach._lease is not None  # noqa: SLF001
        assert len({owner_a._lease.dir, owner_b._lease.dir, attach._lease.dir}) == 3  # noqa: SLF001

        attach_dir = attach._lease.dir  # noqa: SLF001
        owner_a_dir = owner_a._lease.dir  # noqa: SLF001
        owner_b_dir = owner_b._lease.dir  # noqa: SLF001
        attach_channel = cast(Any, attach.channel)

        owner_a.close()
        assert reaped == [owner_a_ep.pid]
        assert not owner_a_dir.exists()
        assert owner_b_dir.exists()
        assert attach_dir.exists()
        assert attach_channel.closed == 0

        owner_b.close()
        assert reaped == [owner_a_ep.pid, owner_b_ep.pid]
        assert not owner_b_dir.exists()
        assert attach_dir.exists()
        assert attach_channel.closed == 0

        attach.close()
        attach_peer.close()
        assert attach_channel.closed == 1
        assert reaped == [owner_a_ep.pid, owner_b_ep.pid]
        assert attach_dir.exists()
    finally:
        FakeChannel.last = None


@pytest.mark.parametrize("action", ["close", "abort"])
def test_detach_launcher_close_and_abort_skip_reap(launcher_env, monkeypatch, action):
    binary, _, endpoint, calls = launcher_env
    spawn_kwargs = {}

    def fake_detach_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs):
        calls["spawn"] += 1
        spawn_kwargs.update(kwargs)
        return FakeProcess()

    monkeypatch.setattr(launcher_module, "spawn_core", fake_detach_spawn)

    launcher = CoreLauncher(binary, detach=True, idle_timeout=2.0)
    assert launcher.start() == endpoint
    assert spawn_kwargs == {"detach": True, "idle_timeout": 2.0}

    channel = cast(Any, launcher.channel)
    getattr(launcher, action)()

    assert channel.closed == 1
    assert calls["reap"] == 0
    assert calls["cleanup"] == 0
    assert launcher.exit_code is None


def test_process_exited_is_false_in_detached_mode(launcher_env, monkeypatch):
    binary, _, _, _ = launcher_env
    proc = FakeProcess()
    monkeypatch.setattr(launcher_module, "spawn_core", lambda *args, **kwargs: proc)

    launcher = CoreLauncher(binary, detach=True, idle_timeout=2.0)
    launcher.start()
    proc.returncode = 1

    assert launcher.process_exited is False


def test_detach_launcher_reaps_process_when_start_fails(launcher_env, monkeypatch):
    binary, _, _, calls = launcher_env

    def fake_detach_spawn(binary_path: Path, runtime: RuntimeLease, **kwargs):
        calls["spawn"] += 1
        return FakeProcess()

    def fail_open(endpoint: Endpoint, timeout: float = 0.0):
        raise CoreGoError(Stage.CHANNEL, ReasonCode.CHANNEL_TIMEOUT, "channel boom")

    monkeypatch.setattr(launcher_module, "spawn_core", fake_detach_spawn)
    monkeypatch.setattr(FakeChannel, "open", fail_open)

    launcher = CoreLauncher(binary, detach=True, idle_timeout=2.0)
    with pytest.raises(CoreGoError) as caught:
        launcher.start()

    assert caught.value.reason is ReasonCode.CHANNEL_TIMEOUT
    assert calls["reap"] == 1
    assert calls["cleanup"] == 1


def test_detach_launcher_rejects_negative_idle_timeout(tmp_path):
    with pytest.raises(ValueError):
        CoreLauncher(tmp_path / "swanlab-core", detach=True, idle_timeout=-0.5)
