"""
@author: Nexisato
@file: base.py
@time: 2026/11/26
@description: CoreGo 与 CoreSyncGo 共享的进程生命周期、状态锁与错误记录基类。
"""

import threading
from enum import Enum
from pathlib import Path
from typing import Any, ClassVar, Generic, Optional, Tuple, TypeVar

from swanlab.sdk.internal.core_go.config import TEARDOWN_TIMEOUT, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc
from swanlab.sdk.internal.pkg import fork, safe

_StateT = TypeVar("_StateT", bound=Enum)


class CoreGoBase(Generic[_StateT]):
    """封装 CoreGo / CoreSyncGo 共用的进程所有权、状态锁、错误归类与传输回收逻辑。"""

    _state_new: ClassVar[Any]
    _state_closed: ClassVar[Any]
    _state_closing: ClassVar[Any]
    _state_aborting: ClassVar[Any]
    _close_abort_states: ClassVar[Tuple[Any, ...]]
    _teardown_log_prefix: ClassVar[str]

    _binary: Optional[Path]
    _detach: bool
    _lock: threading.RLock
    _owner_pid: int
    _state: _StateT
    _launcher: Optional[CoreLauncher]
    _last_error: Optional[CoreGoError]
    _exit_code: Optional[int]

    def _init_lifecycle(
        self,
        *,
        binary: Optional[Path],
        detach: bool,
        initial_state: _StateT,
    ) -> None:
        self._binary = Path(binary) if binary is not None else None
        self._detach = detach
        self._lock = threading.RLock()
        self._owner_pid = fork.current_pid()
        self._state = initial_state
        self._launcher = None
        self._last_error = None
        self._exit_code = None

    @property
    def last_error(self) -> Optional[CoreGoError]:
        return self._last_error

    @property
    def exit_code(self) -> Optional[int]:
        """最近一次 close/abort 回收的 core 进程退出码。"""
        return self._exit_code

    def _teardown_rpc_locked(self) -> Optional[CoreServiceRpc]:
        raise NotImplementedError

    def _clear_transport_locked(self) -> None:
        raise NotImplementedError

    def _record_start_error(self, stage: Stage, reason: ReasonCode, message: str) -> CoreGoError:
        """构造并记录启动失败错误；调用方必须持有 self._lock。"""
        error = CoreGoError(stage, reason, message)
        self._last_error = error
        return error

    def _record_error(self, exc: BaseException, stage: Stage, *, abort: bool = False) -> CoreGoError:
        error = self._as_core_error(exc, stage)
        with self._lock:
            self._last_error = error
        if abort:
            self.abort(error)
        return error

    @staticmethod
    def _as_core_error(exc: BaseException, stage: Stage) -> CoreGoError:
        if isinstance(exc, CoreGoError):
            return exc
        return CoreGoError(stage, ReasonCode.GRPC_ERROR, str(exc), cause=exc)

    def _owned_by_current_process(self) -> bool:
        # 只拦真正的 os.fork()；macOS 默认 spawn，不会继承本对象。
        return not fork.is_forked(self._owner_pid)

    def close(self) -> None:
        """正常关闭当前持有的 Core 传输与子进程；幂等且不抛出。"""
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (self._state_closed, self._state_closing, self._state_aborting):
                return
            if self._state in self._close_abort_states:
                abort_startup = True
                rpc = None
                launcher = None
            elif self._state is self._state_new:
                self._state = self._state_closed
                return
            else:
                abort_startup = False
                self._state = self._state_closing
                rpc = self._teardown_rpc_locked()
                launcher = self._launcher

        if abort_startup:
            self.abort()
            return

        self._release_transport(rpc, launcher, abort_launcher=False)

    def abort(self, cause: Optional[BaseException] = None) -> None:
        """失败或中断路径回收 Core 传输与子进程；幂等且不抛出。"""
        _ = cause
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (self._state_aborting, self._state_closing, self._state_closed):
                return
            self._state = self._state_aborting
            rpc = self._teardown_rpc_locked()
            launcher = self._launcher

        self._release_transport(rpc, launcher, abort_launcher=True)

    def _release_transport(
        self,
        rpc: Optional[CoreServiceRpc],
        launcher: Optional[CoreLauncher],
        *,
        abort_launcher: bool,
    ) -> None:
        """释放 RPC/Launcher 资源并统一收敛到本地 CLOSED；永不抛出。"""
        action = "abort" if abort_launcher else "close"
        try:
            owns_process = getattr(launcher, "owns_process", True) if launcher is not None else True
            if not self._detach and owns_process and rpc is not None:
                with safe.block(level="debug", message=f"{self._teardown_log_prefix} teardown failed during {action}"):
                    rpc.teardown(timeout=TEARDOWN_TIMEOUT)
            if launcher is not None:
                with safe.block(level="debug", message=f"{self._teardown_log_prefix} launcher {action} failed"):
                    if abort_launcher:
                        launcher.abort()
                    else:
                        launcher.close()
                self._exit_code = launcher.exit_code
        finally:
            with self._lock:
                self._clear_transport_locked()
                self._launcher = None
                self._state = self._state_closed
