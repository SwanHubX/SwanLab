"""
@author: Nexisato
@file: launcher.py
@time: 2026/10/27 15:40
@description: CoreLauncher：core 服务进程的本地所有权。

负责编排 runtime →（复用存活实例或持锁按需拉起）→ 等 port-file →
打开 gRPC channel → 幂等 close/abort 回收。
只面向本地端点，不构造业务 RPC，不承担重试策略。
"""

import subprocess
import threading
from enum import Enum
from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.core_go.config import (
    CHANNEL_READY_TIMEOUT,
    DETACH_IDLE_TIMEOUT,
    PORT_FILE_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.process import reap_process, spawn_core, wait_endpoint
from swanlab.sdk.internal.core_go.rpc import Channel
from swanlab.sdk.internal.core_go.runtime import (
    Endpoint,
    RuntimeLease,
    acquire_runtime,
    acquire_start_lock,
    cleanup_runtime,
    discover_endpoint,
    validate_core_binary,
)
from swanlab.sdk.internal.pkg import fork

__all__ = ["CoreLauncher"]


class _LauncherState(str, Enum):
    NEW = "new"
    STARTING = "starting"
    RUNNING = "running"
    CLOSED = "closed"


class CoreLauncher:
    """core 进程生命周期与本地传输通道的管理者。

    ``lease`` 仅用于测试注入（例如预置 core.sock 触发 TCP 回退），
    生产调用不传，由 start() 自行 acquire_runtime()。
    """

    def __init__(
        self,
        binary: Path,
        *,
        lease: Optional[RuntimeLease] = None,
        detach: bool = False,
        idle_timeout: float = DETACH_IDLE_TIMEOUT,
    ) -> None:
        if idle_timeout < 0:
            raise ValueError("idle_timeout must be non-negative")
        self._binary = Path(binary)
        self._injected_lease = lease
        self._detach = detach
        self._idle_timeout = idle_timeout
        self._lease: Optional[RuntimeLease] = None
        self._proc: Optional[subprocess.Popen] = None
        self._channel: Optional[Channel] = None
        self._endpoint: Optional[Endpoint] = None
        self._exit_code: Optional[int] = None
        self._owns_process = False
        self._lock = threading.Lock()
        self._state = _LauncherState.NEW
        self._owner_pid = fork.current_pid()

    @property
    def endpoint(self) -> Endpoint:
        if self._state is _LauncherState.CLOSED:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")
        if self._endpoint is None:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "core launcher has not started")
        return self._endpoint

    @property
    def channel(self) -> Channel:
        if self._state is _LauncherState.CLOSED:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")
        if self._channel is None:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "core launcher has not started")
        return self._channel

    @property
    def exit_code(self) -> Optional[int]:
        """回收后的进程退出码；未启动或尚未回收时为 None。"""
        return self._exit_code

    @property
    def owns_process(self) -> bool:
        """当前 launcher 是否拉起并拥有该 core 子进程；复用已有端点时为 False。"""
        return self._owns_process

    def start(self) -> Endpoint:
        """拉起 core 并等待 channel-ready。

        失败时先回收已创建的 runtime / 进程 / channel，再原样抛出异常；
        成功后再调用 start() 是幂等的，直接返回现有端点。
        """
        if not self._owned_by_current_process():
            raise CoreGoError(
                Stage.LIFECYCLE,
                ReasonCode.FORK_UNSUPPORTED,
                "core launcher belongs to another process",
            )

        with self._lock:
            if self._state is _LauncherState.RUNNING:
                assert self._endpoint is not None
                return self._endpoint
            if self._state is _LauncherState.CLOSED:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")

            self._state = _LauncherState.STARTING
            try:
                validate_core_binary(self._binary)
                lease = (
                    self._injected_lease if self._injected_lease is not None else acquire_runtime(shared=self._detach)
                )
                self._lease = lease
                self._validate_lease(lease)

                # Attach 模式使用每节点每用户共享 runtime，支持多 Rank / 多 Run 复用单个 detached Core；
                # Owner 模式每 Run 独占隔离目录和私有子进程，随 Run teardown 并清理，互不干扰。
                if self._injected_lease is None and self._detach:
                    self._start_shared_locked(lease)
                else:
                    self._spawn_and_connect_locked(lease)
            except BaseException:
                self._cleanup_locked()
                raise

            self._state = _LauncherState.RUNNING
            assert self._endpoint is not None
            return self._endpoint

    def _start_shared_locked(self, lease: RuntimeLease) -> None:
        """先探测并复用已有存活 Core；未命中时持启动锁二次检查再按需拉起。"""
        if self._try_reuse_existing(lease):
            return

        with acquire_start_lock(lease.dir):
            if self._try_reuse_existing(lease):
                return
            self._spawn_and_connect_locked(lease)

    def _try_reuse_existing(self, lease: RuntimeLease) -> bool:
        existing = discover_endpoint(lease.port_file, lease.dir)
        if existing is None:
            return False
        try:
            self._reuse_endpoint_locked(existing)
        except CoreGoError as exc:
            if exc.stage is not Stage.CHANNEL:
                raise
            # 复用可能恰好撞上旧 Core 空闲退出：清空半成品状态并交由后续流程按需重起。
            self._owns_process = False
            self._endpoint = None
            self._channel = None
            return False
        return True

    def _reuse_endpoint_locked(self, endpoint: Endpoint) -> None:
        self._owns_process = False
        self._endpoint = endpoint
        self._channel = Channel.open(endpoint, timeout=CHANNEL_READY_TIMEOUT)

    def _spawn_and_connect_locked(self, lease: RuntimeLease) -> None:
        self._owns_process = True
        if self._injected_lease is None:
            # 持锁且确认旧端点已失效后清理残留 port-file，防止 wait_endpoint 误读旧 PID。
            lease.port_file.unlink(missing_ok=True)
        proc = spawn_core(self._binary, lease, detach=self._detach, idle_timeout=self._idle_timeout)
        self._proc = proc
        try:
            endpoint = wait_endpoint(proc, lease, timeout=PORT_FILE_TIMEOUT)
        except CoreGoError as exc:
            if exc.reason is not ReasonCode.ENDPOINT_BUSY:
                raise
            # 极端并发下被外部启动的实例抢先占住端点，回收失败子进程并尝试复用赢家。
            reap_process(proc)
            self._proc = None
            self._owns_process = False
            if self._try_reuse_existing(lease):
                return
            raise

        self._endpoint = endpoint
        self._channel = Channel.open(endpoint, timeout=CHANNEL_READY_TIMEOUT)

    def close(self) -> None:
        """正常关闭：channel → reap → runtime 清理；不发送任何业务 RPC。"""
        self._shutdown()

    def abort(self) -> None:
        """失败路径回收：与 close 同一套动作，幂等且不抛出。"""
        self._shutdown()

    def _shutdown(self) -> None:
        # fork child 不持有父进程的锁，也不应回收父进程的 core。
        if not self._owned_by_current_process():
            return
        with self._lock:
            if self._state is _LauncherState.CLOSED:
                return
            self._cleanup_locked()

    def _cleanup_locked(self) -> None:
        """调用方必须持有 self._lock；所有回收动作都不允许覆盖原始异常。"""
        startup_failed = self._state is _LauncherState.STARTING
        channel = self._channel
        proc = self._proc
        lease = self._lease
        owns_process = self._owns_process

        self._channel = None
        self._proc = None
        self._lease = None
        self._endpoint = None
        self._owns_process = False
        self._state = _LauncherState.CLOSED

        if channel is not None:
            channel.close()
        # Attach/Detached 模式下，一旦 start() 成功，关闭客户端连接不得杀掉共享 Core；
        # 仅在 start() 本身未就绪前失败时才回收刚拉起的异常子进程。
        reap_owned = owns_process and proc is not None and (not self._detach or startup_failed)
        if reap_owned:
            self._exit_code = reap_process(proc)
        if owns_process and lease is not None and (not self._detach or startup_failed):
            cleanup_runtime(lease)

    def _validate_lease(self, lease: RuntimeLease) -> None:
        if lease.owner_pid != self._owner_pid:
            raise CoreGoError(
                Stage.RUNTIME,
                ReasonCode.RUNTIME_UNAVAILABLE,
                f"runtime lease belongs to pid {lease.owner_pid}, current pid {self._owner_pid}",
            )
        if lease.dir.is_symlink() or not lease.dir.is_dir():
            raise CoreGoError(
                Stage.RUNTIME,
                ReasonCode.RUNTIME_UNAVAILABLE,
                f"runtime dir is missing or not a real directory: {lease.dir}",
            )

    def _owned_by_current_process(self) -> bool:
        # pkg.fork 只做 PID 比较，不依赖 multiprocessing 的 fork/spawn 选择；
        # macOS 默认 spawn，spawn child 不会继承本对象，这里只拦真正的 os.fork()。
        return not fork.is_forked(self._owner_pid)
