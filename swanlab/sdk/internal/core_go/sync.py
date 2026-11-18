"""
@author: Nexisato
@file: sync.py
@time: 2026/11/18
@description: CoreSyncProtocol 的 Go core 适配层，负责独占进程拉起、Sync RPC 转发与生命周期回收。
"""

import threading
from enum import Enum
from pathlib import Path
from typing import Optional

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
    DeliverSyncStartResponse,
)
from swanlab.sdk.internal.core_go.config import (
    TEARDOWN_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc, CoreSyncRpc
from swanlab.sdk.internal.pkg import fork, safe
from swanlab.sdk.protocol.core import CoreSyncProtocol

__all__ = ["CoreSyncGo"]


class _SyncState(str, Enum):
    NEW = "new"
    STARTING = "starting"
    STARTED = "started"
    FLUSHING = "flushing"
    RUNNING = "running"
    FINISHING = "finishing"
    FINISHED = "finished"
    CLOSING = "closing"
    ABORTING = "aborting"
    CLOSED = "closed"


class CoreSyncGo(CoreSyncProtocol):
    """CoreSyncProtocol 的 Go core 实现，以 Owner Mode 独占一个 swanlab-core 进程完成单次 sync。

    状态机：NEW → STARTING → STARTED → FLUSHING → RUNNING → FINISHING → FINISHED → CLOSED；
    启动/Flush 失败或外部中断走 ABORTING → CLOSED。
    """

    def __init__(self, *, binary: Optional[Path] = None) -> None:
        super().__init__()
        self._binary = Path(binary) if binary is not None else None
        self._lock = threading.RLock()
        self._owner_pid = fork.current_pid()
        self._state = _SyncState.NEW
        self._launcher: Optional[CoreLauncher] = None
        self._core_rpc: Optional[CoreServiceRpc] = None
        self._sync_rpc: Optional[CoreSyncRpc] = None
        self._last_error: Optional[CoreGoError] = None
        self._abort_cause: Optional[BaseException] = None
        self._exit_code: Optional[int] = None

    @property
    def last_error(self) -> Optional[CoreGoError]:
        return self._last_error

    @property
    def exit_code(self) -> Optional[int]:
        """最近一次 close/abort 回收的 core 进程退出码。"""
        return self._exit_code

    # ---------------------------------- 启动与 Flush ----------------------------------

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        with self._lock:
            if self._state is _SyncState.CLOSED:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core is closed")
            if self._state is not _SyncState.NEW:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.ALREADY_STARTED, "sync has already been started")
            binary = self._binary
            if binary is None:
                error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "core binary is not configured")
                self._last_error = error
                self._state = _SyncState.CLOSED
                return DeliverSyncStartResponse(success=False, message=str(error))
            self._state = _SyncState.STARTING

        try:
            sync_rpc = self._connect_and_spinup(binary)
            response = sync_rpc.deliver_sync_start(start_request)
            if not response.success:
                raise CoreGoError(
                    Stage.RUN_START,
                    ReasonCode.RUN_START_REJECTED,
                    response.message or "DeliverSyncStart rejected",
                )

            with self._lock:
                if self._state is not _SyncState.STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
                self._state = _SyncState.STARTED
            return response
        except BaseException as exc:
            error = self._record_error(exc, Stage.LIFECYCLE, abort=True)
            if not isinstance(exc, Exception):
                raise
            return DeliverSyncStartResponse(success=False, message=str(error))

    def _connect_and_spinup(self, binary: Path) -> CoreSyncRpc:
        launcher = CoreLauncher(binary, detach=False)
        with self._lock:
            if self._state is not _SyncState.STARTING:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
            self._launcher = launcher
        launcher.start()

        core_rpc = CoreServiceRpc(launcher.channel)
        sync_rpc = CoreSyncRpc(launcher.channel)
        with self._lock:
            if self._state is not _SyncState.STARTING:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
            self._core_rpc = core_rpc
            self._sync_rpc = sync_rpc

        core_rpc.spinup()
        return sync_rpc

    def _fail_start(self, stage: Stage, reason: ReasonCode, message: str) -> DeliverSyncStartResponse:
        """调用方必须持有 self._lock。"""
        error = CoreGoError(stage, reason, message)
        self._last_error = error
        return DeliverSyncStartResponse(success=False, message=str(error))

    def deliver_sync_flush(self) -> DeliverSyncFlushResponse:
        with self._lock:
            if self._state is not _SyncState.STARTED or self._sync_rpc is None:
                if self._state is _SyncState.CLOSED:
                    reason = ReasonCode.CLOSED
                elif self._state in (
                    _SyncState.FLUSHING,
                    _SyncState.RUNNING,
                    _SyncState.FINISHING,
                    _SyncState.FINISHED,
                ):
                    reason = ReasonCode.ALREADY_STARTED
                else:
                    reason = ReasonCode.NOT_STARTED
                error = CoreGoError(Stage.LIFECYCLE, reason, "sync is not in started state; refusing to flush")
                self._last_error = error
                return DeliverSyncFlushResponse(success=False, message=str(error))
            sync_rpc = self._sync_rpc
            self._state = _SyncState.FLUSHING

        try:
            response = sync_rpc.deliver_sync_flush()
            if not response.success:
                raise CoreGoError(
                    Stage.RUN_START,
                    ReasonCode.RUN_START_REJECTED,
                    response.message or "DeliverSyncFlush rejected",
                )
            with self._lock:
                if self._state is not _SyncState.FLUSHING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during flush")
                self._state = _SyncState.RUNNING
            return response
        except BaseException as exc:
            error = self._record_error(exc, Stage.RPC, abort=True)
            if not isinstance(exc, Exception):
                raise
            return DeliverSyncFlushResponse(success=False, message=str(error))

    # ---------------------------------- 进度查询与收尾 ----------------------------------

    def get_operation_stats(self) -> GetOperationStatsResponse:
        with self._lock:
            sync_rpc = self._sync_rpc
            state = self._state
        if sync_rpc is None or state not in (_SyncState.RUNNING, _SyncState.FINISHING, _SyncState.FINISHED):
            reason = ReasonCode.CLOSED if state is _SyncState.CLOSED else ReasonCode.NOT_STARTED
            error = CoreGoError(Stage.LIFECYCLE, reason, "sync session is not running")
            with self._lock:
                self._last_error = error
            return GetOperationStatsResponse(success=False, message=str(error))

        try:
            return sync_rpc.get_operation_stats()
        except Exception as exc:
            error = self._record_error(exc, Stage.RPC)
            return GetOperationStatsResponse(success=False, message=str(error))

    def confirm_sync_finish(self) -> ConfirmSyncFinishResponse:
        with self._lock:
            if self._state is not _SyncState.RUNNING or self._sync_rpc is None:
                reason = ReasonCode.CLOSED if self._state is _SyncState.CLOSED else ReasonCode.NOT_STARTED
                error = CoreGoError(Stage.LIFECYCLE, reason, "sync is not running; refusing to confirm")
                self._last_error = error
                return ConfirmSyncFinishResponse(success=False, message=str(error))
            sync_rpc = self._sync_rpc
            self._state = _SyncState.FINISHING

        try:
            response = sync_rpc.confirm_sync_finish()
        except BaseException as exc:
            if not isinstance(exc, Exception):
                self._record_error(exc, Stage.RPC, abort=True)
                raise
            error = self._record_error(exc, Stage.RPC)
            with self._lock:
                if self._state is _SyncState.FINISHING:
                    self._state = _SyncState.FINISHED
            return ConfirmSyncFinishResponse(success=False, message=str(error))

        with self._lock:
            if self._state is _SyncState.FINISHING:
                self._state = _SyncState.FINISHED
            if not response.success:
                self._last_error = CoreGoError(
                    Stage.LIFECYCLE,
                    ReasonCode.RUN_FINISH_REJECTED,
                    response.message or "ConfirmSyncFinish rejected",
                )
        return response

    # ---------------------------------- 生命周期回收 ----------------------------------

    def close(self) -> None:
        """正常关闭：Owner Mode 下 teardown + 回收独占进程；幂等且不抛出。"""
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (_SyncState.CLOSED, _SyncState.CLOSING, _SyncState.ABORTING):
                return
            if self._state in (_SyncState.STARTING, _SyncState.FLUSHING):
                abort_startup = True
                core_rpc = None
                launcher = None
            elif self._state is _SyncState.NEW:
                self._state = _SyncState.CLOSED
                return
            else:
                abort_startup = False
                self._state = _SyncState.CLOSING
                core_rpc = self._core_rpc
                launcher = self._launcher

        if abort_startup:
            self.abort()
            return

        self._release_transport(core_rpc, launcher, abort_launcher=False)

    def abort(self, cause: Optional[BaseException] = None) -> None:
        """失败或中断路径回收：Owner Mode 下尽力 teardown + 有界回收独占进程；幂等且不抛出。"""
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (_SyncState.ABORTING, _SyncState.CLOSING, _SyncState.CLOSED):
                return
            self._state = _SyncState.ABORTING
            if cause is not None:
                self._abort_cause = cause
            core_rpc = self._core_rpc
            launcher = self._launcher

        self._release_transport(core_rpc, launcher, abort_launcher=True)

    def _release_transport(
        self,
        core_rpc: Optional[CoreServiceRpc],
        launcher: Optional[CoreLauncher],
        *,
        abort_launcher: bool,
    ) -> None:
        """释放 RPC/Launcher 资源并统一收敛到本地 CLOSED；永不抛出。"""
        action = "abort" if abort_launcher else "close"
        try:
            owns_process = getattr(launcher, "owns_process", True) if launcher is not None else True
            if owns_process and core_rpc is not None:
                with safe.block(level="debug", message=f"sync core teardown failed during {action}"):
                    core_rpc.teardown(timeout=TEARDOWN_TIMEOUT)
            if launcher is not None:
                with safe.block(level="debug", message=f"sync core launcher {action} failed"):
                    if abort_launcher:
                        launcher.abort()
                    else:
                        launcher.close()
                self._exit_code = launcher.exit_code
        finally:
            with self._lock:
                self._core_rpc = None
                self._sync_rpc = None
                self._launcher = None
                self._state = _SyncState.CLOSED

    def _owned_by_current_process(self) -> bool:
        return not fork.is_forked(self._owner_pid)

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
