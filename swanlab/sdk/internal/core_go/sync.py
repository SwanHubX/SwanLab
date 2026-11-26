"""
@author: Nexisato
@file: sync.py
@time: 2026/11/18
@description: CoreSyncProtocol 的 Go core 适配层，负责独占进程拉起、Sync RPC 转发与生命周期回收。
"""

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
from swanlab.sdk.internal.core_go.base import CoreGoBase
from swanlab.sdk.internal.core_go.config import (
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc, CoreSyncRpc
from swanlab.sdk.protocol.core import CoreSyncProtocol

__all__ = ["CoreSyncGo"]


class SyncState(str, Enum):
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


class CoreSyncGo(CoreGoBase[SyncState], CoreSyncProtocol):
    """CoreSyncProtocol 的 Go core 实现，以 Owner Mode 独占一个 swanlab-core 进程完成单次 sync。

    状态机：NEW → STARTING → STARTED → FLUSHING → RUNNING → FINISHING → FINISHED → CLOSED；
    启动/Flush 失败或外部中断走 ABORTING → CLOSED。
    """

    _state_new = SyncState.NEW
    _state_closed = SyncState.CLOSED
    _state_closing = SyncState.CLOSING
    _state_aborting = SyncState.ABORTING
    _close_abort_states = (SyncState.STARTING, SyncState.FLUSHING)
    _teardown_log_prefix = "sync core"

    def __init__(self, *, binary: Optional[Path] = None) -> None:
        CoreSyncProtocol.__init__(self)
        self._init_lifecycle(binary=binary, detach=False, initial_state=SyncState.NEW)
        self._core_rpc: Optional[CoreServiceRpc] = None
        self._sync_rpc: Optional[CoreSyncRpc] = None

    def _teardown_rpc_locked(self) -> Optional[CoreServiceRpc]:
        return self._core_rpc

    def _clear_transport_locked(self) -> None:
        self._core_rpc = None
        self._sync_rpc = None

    # ---------------------------------- 启动与 Flush ----------------------------------

    def deliver_sync_start(self, start_request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        with self._lock:
            if self._state is SyncState.CLOSED:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core is closed")
            if self._state is not SyncState.NEW:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.ALREADY_STARTED, "sync has already been started")
            binary = self._binary
            if binary is None:
                resp = self._fail_start(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "core binary is not configured")
                self._state = SyncState.CLOSED
                return resp
            self._state = SyncState.STARTING

        try:
            sync_rpc = self._connect_and_spinup(binary)
            response = sync_rpc.deliver_sync_start(start_request)

            with self._lock:
                if self._state is not SyncState.STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
                self._state = SyncState.STARTED
            return response
        except BaseException as exc:
            error = self._record_error(exc, Stage.LIFECYCLE, abort=True)
            if not isinstance(exc, Exception):
                raise
            return DeliverSyncStartResponse(success=False, message=str(error))

    def _connect_and_spinup(self, binary: Path) -> CoreSyncRpc:
        launcher = CoreLauncher(binary, detach=False)
        with self._lock:
            if self._state is not SyncState.STARTING:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
            self._launcher = launcher
        launcher.start()

        core_rpc = CoreServiceRpc(launcher.channel)
        sync_rpc = CoreSyncRpc(launcher.channel)
        with self._lock:
            if self._state is not SyncState.STARTING:
                raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during startup")
            self._core_rpc = core_rpc
            self._sync_rpc = sync_rpc

        core_rpc.spinup()
        return sync_rpc

    def _fail_start(self, stage: Stage, reason: ReasonCode, message: str) -> DeliverSyncStartResponse:
        """调用方必须持有 self._lock。"""
        error = self._record_start_error(stage, reason, message)
        return DeliverSyncStartResponse(success=False, message=str(error))

    def deliver_sync_flush(self) -> DeliverSyncFlushResponse:
        with self._lock:
            if self._state is not SyncState.STARTED or self._sync_rpc is None:
                if self._state is SyncState.CLOSED:
                    reason = ReasonCode.CLOSED
                elif self._state in (
                    SyncState.FLUSHING,
                    SyncState.RUNNING,
                    SyncState.FINISHING,
                    SyncState.FINISHED,
                ):
                    reason = ReasonCode.ALREADY_STARTED
                else:
                    reason = ReasonCode.NOT_STARTED
                error = CoreGoError(Stage.LIFECYCLE, reason, "sync is not in started state; refusing to flush")
                self._last_error = error
                return DeliverSyncFlushResponse(success=False, message=str(error))
            sync_rpc = self._sync_rpc
            self._state = SyncState.FLUSHING

        try:
            response = sync_rpc.deliver_sync_flush()
            with self._lock:
                if self._state is not SyncState.FLUSHING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "sync core was closed during flush")
                self._state = SyncState.RUNNING
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
        if sync_rpc is None or state not in (SyncState.RUNNING, SyncState.FINISHING, SyncState.FINISHED):
            reason = ReasonCode.CLOSED if state is SyncState.CLOSED else ReasonCode.NOT_STARTED
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
            if self._state is not SyncState.RUNNING or self._sync_rpc is None:
                reason = ReasonCode.CLOSED if self._state is SyncState.CLOSED else ReasonCode.NOT_STARTED
                error = CoreGoError(Stage.LIFECYCLE, reason, "sync is not running; refusing to confirm")
                self._last_error = error
                return ConfirmSyncFinishResponse(success=False, message=str(error))
            sync_rpc = self._sync_rpc
            self._state = SyncState.FINISHING

        try:
            response = sync_rpc.confirm_sync_finish()
        except BaseException as exc:
            if not isinstance(exc, Exception):
                self._record_error(exc, Stage.RPC, abort=True)
                raise
            error = self._record_error(exc, Stage.RPC)
            with self._lock:
                if self._state is SyncState.FINISHING:
                    self._state = SyncState.FINISHED
            return ConfirmSyncFinishResponse(success=False, message=str(error))

        with self._lock:
            if self._state is SyncState.FINISHING:
                self._state = SyncState.FINISHED
        return response
