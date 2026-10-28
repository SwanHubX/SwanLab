"""
@author: Nexisato
@file: launcher.py
@time: 2026/10/27 15:40
@description: CoreLauncher：core 服务进程的本地所有权。

唯一同时持有 process、runtime 与 channel 的对象，协议无关：
start() 负责 runtime → spawn → 等待端点 → channel-ready；
close() / abort() 负责分级回收并删除 runtime，Teardown 由 CoreGo 在调用前完成。
不认识 CoreService，也不发送任何业务 RPC；进程机制本身在 process/ 下。
超时与文件名取自 config.py，不在本模块散落硬编码值。
"""

import subprocess
import threading
from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.core_go.config import CHANNEL_READY_TIMEOUT, PORT_FILE_TIMEOUT
from swanlab.sdk.internal.core_go.discovery import Endpoint
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.process import (
    RuntimeLease,
    acquire_runtime,
    cleanup_runtime,
    reap_process,
    spawn_core,
    validate_core_binary,
    wait_endpoint,
)
from swanlab.sdk.internal.core_go.rpc import Channel
from swanlab.sdk.internal.pkg import fork

__all__ = ["CoreLauncher"]

_NEW = "new"
_STARTING = "starting"
_RUNNING = "running"
_CLOSED = "closed"


class CoreLauncher:
    """一次 core 进程生命周期的唯一所有者。

    ``lease`` 仅用于测试注入（例如预置 core.sock 触发 TCP 回退），
    生产调用不传，由 start() 自行 acquire_runtime()。
    """

    def __init__(self, binary: Path, *, lease: Optional[RuntimeLease] = None) -> None:
        self._binary = Path(binary)
        self._injected_lease = lease
        self._lease: Optional[RuntimeLease] = None
        self._proc: Optional[subprocess.Popen] = None
        self._channel: Optional[Channel] = None
        self._endpoint: Optional[Endpoint] = None
        self._lock = threading.Lock()
        self._state = _NEW
        self._owner_pid = fork.current_pid()

    @property
    def endpoint(self) -> Endpoint:
        if self._state == _CLOSED:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")
        if self._endpoint is None:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "core launcher has not started")
        return self._endpoint

    @property
    def channel(self) -> Channel:
        if self._state == _CLOSED:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")
        if self._channel is None:
            raise CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "core launcher has not started")
        return self._channel

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
            if self._state == _RUNNING:
                assert self._endpoint is not None
                return self._endpoint
            if self._state == _CLOSED:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core launcher is closed")

            self._state = _STARTING
            try:
                validate_core_binary(self._binary)
                lease = self._injected_lease if self._injected_lease is not None else acquire_runtime()
                self._lease = lease
                self._validate_lease(lease)
                self._proc = spawn_core(self._binary, lease)
                self._endpoint = wait_endpoint(self._proc, lease, timeout=PORT_FILE_TIMEOUT)
                self._channel = Channel.open(self._endpoint, timeout=CHANNEL_READY_TIMEOUT)
            except BaseException:
                self._cleanup_locked()
                raise

            self._state = _RUNNING
            return self._endpoint

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
            if self._state == _CLOSED:
                return
            self._cleanup_locked()

    def _cleanup_locked(self) -> None:
        """调用方必须持有 self._lock；所有回收动作都不允许覆盖原始异常。"""
        channel = self._channel
        proc = self._proc
        lease = self._lease

        self._channel = None
        self._proc = None
        self._lease = None
        self._endpoint = None
        self._state = _CLOSED

        if channel is not None:
            channel.close()
        if proc is not None:
            reap_process(proc)
        if lease is not None:
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
