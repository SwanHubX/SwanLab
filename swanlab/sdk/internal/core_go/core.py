"""
@author: Nexisato
@file: core.py
@time: 2026/10/27 15:40
@description: CoreProtocol 的 Go core 适配层，负责进程拉起、RPC 转发与生命周期回收。
"""

from enum import Enum
from pathlib import Path
from typing import Callable, List, Optional, cast

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord, MediaRecord, ScalarRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.core_go.base import _CoreGoBase
from swanlab.sdk.internal.core_go.config import (
    DETACH_IDLE_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
    classify_start_error,
)
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc
from swanlab.sdk.protocol import CoreProtocol, NullCoreProtocol
from swanlab.sdk.typings.run import ModeType

__all__ = ["CoreGo"]


class _CoreState(str, Enum):
    NEW = "new"
    STARTING = "starting"
    RUNNING = "running"
    FINISHED = "finished"
    CLOSING = "closing"
    ABORTING = "aborting"
    CLOSED = "closed"


_FORK_MESSAGE = "CoreGo is detached in a forked process; run RPCs are unsupported"


class _DetachedCoreGo(NullCoreProtocol):
    """fork child 的本地占位：不持有父进程资源，也不 reconnect / attach / 回收父 Core。"""

    def __init__(self, mode: ModeType, binary: Optional[Path]) -> None:
        super().__init__(mode, message=_FORK_MESSAGE)
        self._binary = binary
        self._last_error = CoreGoError(Stage.FORK, ReasonCode.FORK_UNSUPPORTED, _FORK_MESSAGE)

    @property
    def last_error(self) -> CoreGoError:
        return self._last_error

    @staticmethod
    def _start_when_disabled(start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return DeliverRunStartResponse(success=False, message=_FORK_MESSAGE)


class CoreGo(_CoreGoBase[_CoreState], CoreProtocol):
    """CoreProtocol 的 Go core 实现，一次实例只绑定一个 run session。

    状态机：NEW → STARTING → RUNNING → FINISHED → CLOSED；失败走 ABORTING。
    """

    _state_new = _CoreState.NEW
    _state_closed = _CoreState.CLOSED
    _state_closing = _CoreState.CLOSING
    _state_aborting = _CoreState.ABORTING
    _close_abort_states = (_CoreState.STARTING,)
    _teardown_log_prefix = "core"

    # upsert 失败交给 consumer 分类重试或终结，不在本层吞掉。
    raises_on_upsert_error = True
    # deliver 失败时会话保留在 core 内等待重投，跳过 confirm。
    confirm_after_deliver_failure = False
    reports_runtime_stats = True

    def __init__(
        self,
        mode: ModeType,
        *,
        binary: Optional[Path] = None,
        detach: bool = False,
        idle_timeout: float = DETACH_IDLE_TIMEOUT,
    ) -> None:
        CoreProtocol.__init__(self, mode)
        if idle_timeout < 0:
            raise ValueError("idle_timeout must be non-negative")
        self._init_lifecycle(binary=binary, detach=detach, initial_state=_CoreState.NEW)
        self._idle_timeout = idle_timeout
        self._rpc: Optional[CoreServiceRpc] = None
        self._session_id: Optional[str] = None

    def _teardown_rpc_locked(self) -> Optional[CoreServiceRpc]:
        return self._rpc

    def _clear_transport_locked(self) -> None:
        self._rpc = None
        self._session_id = None

    # ---------------------------------- 实验开始 ----------------------------------

    def _start_when_enabled(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        with self._lock:
            if self._state is _CoreState.CLOSED:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.CLOSED, "core is closed")
            if self._state is not _CoreState.NEW:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.ALREADY_STARTED, "run has already been started")
            binary = self._binary
            if binary is None:
                resp = self._fail_start(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "core binary is not configured")
                self._state = _CoreState.CLOSED
                return resp
            self._state = _CoreState.STARTING

        try:
            rpc = self._connect_and_spinup(binary)
            response = rpc.deliver_run_start(start_request)

            with self._lock:
                if self._state is not _CoreState.STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._session_id = response.session_id
                self._state = _CoreState.RUNNING
            return response
        except BaseException as exc:
            error = self._record_error(exc, Stage.LIFECYCLE, abort=True)
            # 清理后原样抛出，避免协议层吞掉 KeyboardInterrupt/SystemExit。
            if not isinstance(exc, Exception):
                raise
            return DeliverRunStartResponse(
                success=False,
                message=str(error),
                error_kind=classify_start_error(error),
            )

    def _connect_and_spinup(self, binary: Path) -> CoreServiceRpc:
        """建立通道并完成 Spinup；复用即将空闲退出的 detached Core 时重试拉起。"""
        max_attempts = 3 if self._detach else 1
        for attempt in range(max_attempts):
            launcher = CoreLauncher(binary, detach=self._detach, idle_timeout=self._idle_timeout)
            with self._lock:
                if self._state is not _CoreState.STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._launcher = launcher
            launcher.start()

            rpc = CoreServiceRpc(launcher.channel)
            with self._lock:
                if self._state is not _CoreState.STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._rpc = rpc

            try:
                rpc.spinup()
                return rpc
            except CoreGoError as exc:
                retryable_shutdown = (
                    self._detach
                    and not launcher.owns_process
                    and exc.stage is Stage.SPINUP
                    and exc.reason in (ReasonCode.FAILED_PRECONDITION, ReasonCode.GRPC_UNAVAILABLE)
                )
                if not retryable_shutdown or attempt == max_attempts - 1:
                    raise
                with self._lock:
                    if self._state is _CoreState.STARTING:
                        self._rpc = None
                        self._launcher = None
                launcher.close()
        raise AssertionError("unreachable")

    def _fail_start(self, stage: Stage, reason: ReasonCode, message: str) -> DeliverRunStartResponse:
        """调用方必须持有 self._lock。"""
        error = self._record_start_error(stage, reason, message)
        return DeliverRunStartResponse(
            success=False,
            message=str(error),
            error_kind=classify_start_error(error),
        )

    # ---------------------------------- 数据上报 ----------------------------------

    def _upsert_columns_when_enabled(self, columns: List[ColumnRecord]) -> None:
        self._forward_upsert("columns", lambda rpc, session_id: rpc.upsert_columns(columns, session_id))

    def _upsert_scalars_when_enabled(self, scalars: List[ScalarRecord]) -> None:
        self._forward_upsert("scalars", lambda rpc, session_id: rpc.upsert_scalars(scalars, session_id))

    def _upsert_media_when_enabled(self, media: List[MediaRecord]) -> None:
        self._forward_upsert("media", lambda rpc, session_id: rpc.upsert_media(media, session_id))

    # 直方图与标量走同一 upsert 路由，失败同样由 consumer 分类重试或终结。
    def _upsert_histograms_when_enabled(self, histograms: List[HistogramRecord]) -> None:
        self._forward_upsert("histograms", lambda rpc, session_id: rpc.upsert_histograms(histograms, session_id))

    def _upsert_logs_when_enabled(self, logs: List[LogRecord]) -> None:
        self._forward_upsert("logs", lambda rpc, session_id: rpc.upsert_logs(logs, session_id))

    def _upsert_saves_when_enabled(self, saves: List[SaveRecord]) -> None:
        self._forward_upsert("saves", lambda rpc, session_id: rpc.upsert_saves(saves, session_id))

    def _forward_upsert(self, name: str, call: Callable[[CoreServiceRpc, str], None]) -> None:
        """转发一次 upsert；失败抛 CoreGoError，由 consumer 按 retryable 分类。"""
        rpc: Optional[CoreServiceRpc] = None
        session_id: Optional[str] = None
        with self._lock:
            if self._state is _CoreState.CLOSED:
                error: Optional[CoreGoError] = CoreGoError(
                    Stage.LIFECYCLE, ReasonCode.CLOSED, f"core is closed; refusing to upsert {name}"
                )
            elif self._state is not _CoreState.RUNNING or self._session_id is None:
                error = CoreGoError(
                    Stage.LIFECYCLE, ReasonCode.NOT_STARTED, f"run is not running; refusing to upsert {name}"
                )
            else:
                error = None
                rpc = self._rpc
                session_id = self._session_id
            if error is not None:
                self._last_error = error

        if error is not None:
            raise error

        assert rpc is not None
        assert session_id is not None
        try:
            call(rpc, session_id)
        except Exception as exc:
            failure = self._as_core_error(exc, Stage.RPC)
            if failure.retryable and self._owner_process_exited():
                failure = CoreGoError(
                    failure.stage,
                    failure.reason,
                    f"{failure.message}; owned core process exited",
                    grpc_name=failure.grpc_name,
                    grpc_code=failure.grpc_code,
                    detail=failure.detail,
                    cause=failure,
                    retryable=False,
                )
            with self._lock:
                self._last_error = failure
            raise failure from exc

    def _owner_process_exited(self) -> bool:
        with self._lock:
            launcher = self._launcher
        return launcher is not None and launcher.process_exited

    # ---------------------------------- 通知运行结束 ----------------------------------

    def _finish_when_enabled(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        with self._lock:
            if self._state is not _CoreState.RUNNING or self._session_id is None:
                error = CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "run is not running; refusing to finish")
                self._last_error = error
                return DeliverRunFinishResponse(success=False, message=str(error))
            rpc = self._rpc
            session_id = self._session_id

        assert rpc is not None
        try:
            response = rpc.deliver_run_finish(finish_request, session_id)
        except Exception as exc:
            error = self._as_core_error(exc, Stage.RPC)
            with self._lock:
                self._last_error = error
            return DeliverRunFinishResponse(success=False, message=str(error))

        with self._lock:
            if self._state is _CoreState.RUNNING and self._session_id == session_id:
                self._state = _CoreState.FINISHED
        return response

    # ---------------------------------- 进度查询 ----------------------------------

    def _get_operation_when_online(self) -> GetOperationStatsResponse:
        with self._lock:
            session_id = self._session_id
            rpc = self._rpc
            state = self._state
        if session_id is None or rpc is None:
            reason = ReasonCode.CLOSED if state is _CoreState.CLOSED else ReasonCode.NOT_STARTED
            error = CoreGoError(Stage.LIFECYCLE, reason, "run session is not bound")
            with self._lock:
                self._last_error = error
            return GetOperationStatsResponse(success=False, message=str(error))

        try:
            return rpc.get_operation_stats(session_id)
        except Exception as exc:
            error = self._as_core_error(exc, Stage.RPC)
            with self._lock:
                self._last_error = error
            return GetOperationStatsResponse(success=False, message=str(error))

    # ---------------------------------- 确认完成 ----------------------------------

    def _confirm_finish_when_enabled(self) -> ConfirmRunFinishResponse:
        with self._lock:
            self._last_error = None
            if self._state is not _CoreState.FINISHED or self._session_id is None:
                reason = ReasonCode.CLOSED if self._state is _CoreState.CLOSED else ReasonCode.NOT_STARTED
                error = CoreGoError(Stage.LIFECYCLE, reason, "run is not finished; refusing to confirm")
                self._last_error = error
                return ConfirmRunFinishResponse(success=False, message=str(error))
            rpc = self._rpc
            session_id = self._session_id

        assert rpc is not None
        try:
            response = rpc.confirm_run_finish(session_id)
        except Exception as exc:
            error = self._as_core_error(exc, Stage.RPC)
            with self._lock:
                self._last_error = error
            return ConfirmRunFinishResponse(success=False, message=str(error))

        with self._lock:
            if self._session_id == session_id:
                self._session_id = None
        return response

    # ---------------------------------- 进程 fork ----------------------------------

    def fork(self) -> "CoreProtocol":
        # 不取锁：fork 时锁可能由子进程不存在的线程持有。
        return _DetachedCoreGo(cast(ModeType, self._mode), self._binary)
