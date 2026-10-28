"""
@author: Nexisato
@file: core.py
@time: 2026/10/27 15:40
@description: CoreGo：CoreProtocol 在 Go core 上的适配层。

只做 run 级编排：spinup + DeliverRunStart 保存 session_id，Upsert 薄转发，
finish/stats/confirm 一一映射；本地存储、批量与重试由 Go core 或上游负责。
CoreEnum.CORE（"Core"）对应本模块，CoreEnum.CORE_PYTHON 对应 CorePython。
"""

import threading
from pathlib import Path
from typing import Callable, List, Optional, Set, cast

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.core_go.config import TEARDOWN_TIMEOUT
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc.core_service import CoreServiceRpc
from swanlab.sdk.internal.pkg import console, fork, safe
from swanlab.sdk.protocol import CoreProtocol, NullCoreProtocol
from swanlab.sdk.typings.run import ModeType

__all__ = ["CoreGo"]

_NEW = "new"
_STARTING = "starting"
_RUNNING = "running"
_FINISHED = "finished"
_CLOSING = "closing"
_ABORTING = "aborting"
_CLOSED = "closed"

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


class CoreGo(CoreProtocol):
    """CoreProtocol 的 Go core 实现，一次实例只绑定一个 run session。

    状态机：NEW → STARTING → RUNNING → FINISHED → CLOSED；失败走 ABORTING。
    阻塞的 spawn/RPC/reap 不持锁，close/abort 先切终态再清理。
    """

    def __init__(self, mode: ModeType, *, binary: Optional[Path] = None) -> None:
        super().__init__(mode)
        self._binary = Path(binary) if binary is not None else None
        self._lock = threading.RLock()
        self._owner_pid = fork.current_pid()
        self._state = _NEW
        self._launcher: Optional[CoreLauncher] = None
        self._rpc: Optional[CoreServiceRpc] = None
        self._session_id: Optional[str] = None
        self._last_error: Optional[CoreGoError] = None
        self._abort_cause: Optional[BaseException] = None
        self._warned: Set[str] = set()

    @property
    def last_error(self) -> Optional[CoreGoError]:
        return self._last_error

    # ---------------------------------- 实验开始 ----------------------------------

    def _start_when_online(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        with self._lock:
            if self._state == _CLOSED:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.CLOSED, "core is closed")
            if self._state != _NEW:
                return self._fail_start(Stage.LIFECYCLE, ReasonCode.ALREADY_STARTED, "run has already been started")
            binary = self._binary
            if binary is None:
                error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "core binary is not configured")
                self._last_error = error
                self._state = _CLOSED
                return DeliverRunStartResponse(success=False, message=str(error))
            self._state = _STARTING

        try:
            launcher = CoreLauncher(binary)
            with self._lock:
                if self._state != _STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._launcher = launcher
            launcher.start()

            rpc = CoreServiceRpc(launcher.channel)
            with self._lock:
                if self._state != _STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._rpc = rpc

            rpc.spinup()
            response = rpc.deliver_run_start(start_request)
            if not response.success:
                raise CoreGoError(
                    Stage.RUN_START,
                    ReasonCode.RUN_START_REJECTED,
                    response.message or "DeliverRunStart rejected",
                )
            if not response.session_id:
                raise CoreGoError(
                    Stage.RUN_START,
                    ReasonCode.PROTOCOL_VIOLATION,
                    "DeliverRunStart returned an empty session_id",
                )

            with self._lock:
                if self._state != _STARTING:
                    raise CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core was closed during startup")
                self._session_id = response.session_id
                self._state = _RUNNING
            return response
        except BaseException as exc:
            error = self._as_core_error(exc, Stage.LIFECYCLE)
            with self._lock:
                self._last_error = error
            self.abort(error)
            # 清理后原样抛出，避免协议层吞掉 KeyboardInterrupt/SystemExit。
            if not isinstance(exc, Exception):
                raise
            return DeliverRunStartResponse(success=False, message=str(error))

    def _start_when_local(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._reject_mode("local")

    def _start_when_offline(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._reject_mode("offline")

    def _reject_mode(self, mode: str) -> DeliverRunStartResponse:
        # 模式本身有效，只是当前 Go build 不承载 local/offline。
        error = CoreGoError(
            Stage.LIFECYCLE,
            ReasonCode.MODE_UNSUPPORTED,
            f"core_go only implements online mode in this build; mode={mode!r} is not handled by the Go core",
        )
        with self._lock:
            self._last_error = error
        return DeliverRunStartResponse(success=False, message=str(error))

    def _fail_start(self, stage: Stage, reason: ReasonCode, message: str) -> DeliverRunStartResponse:
        """调用方必须持有 self._lock。"""
        error = CoreGoError(stage, reason, message)
        self._last_error = error
        return DeliverRunStartResponse(success=False, message=str(error))

    # ---------------------------------- 数据上报 ----------------------------------

    def _upsert_columns_when_local(self, columns: List[ColumnRecord]) -> None:
        return None

    def _upsert_columns_when_offline(self, columns: List[ColumnRecord]) -> None:
        return None

    def _upsert_columns_when_online(self, columns: List[ColumnRecord]) -> None:
        self._forward_upsert("columns", lambda rpc, session_id: rpc.upsert_columns(columns, session_id))

    def _upsert_scalars_when_local(self, scalars: List[ScalarRecord]) -> None:
        return None

    def _upsert_scalars_when_offline(self, scalars: List[ScalarRecord]) -> None:
        return None

    def _upsert_scalars_when_online(self, scalars: List[ScalarRecord]) -> None:
        self._forward_upsert("scalars", lambda rpc, session_id: rpc.upsert_scalars(scalars, session_id))

    def _upsert_media_when_local(self, media: List[MediaRecord]) -> None:
        return None

    def _upsert_media_when_offline(self, media: List[MediaRecord]) -> None:
        return None

    def _upsert_media_when_online(self, media: List[MediaRecord]) -> None:
        self._forward_upsert("media", lambda rpc, session_id: rpc.upsert_media(media, session_id))

    def _upsert_logs_when_local(self, logs: List[LogRecord]) -> None:
        return None

    def _upsert_logs_when_offline(self, logs: List[LogRecord]) -> None:
        return None

    def _upsert_logs_when_online(self, logs: List[LogRecord]) -> None:
        self._forward_upsert("logs", lambda rpc, session_id: rpc.upsert_logs(logs, session_id))

    def _upsert_saves_when_local(self, saves: List[SaveRecord]) -> None:
        return None

    def _upsert_saves_when_offline(self, saves: List[SaveRecord]) -> None:
        return None

    def _upsert_saves_when_online(self, saves: List[SaveRecord]) -> None:
        self._forward_upsert("saves", lambda rpc, session_id: rpc.upsert_saves(saves, session_id))

    def _forward_upsert(self, name: str, call: Callable[[CoreServiceRpc, str], None]) -> None:
        rpc: Optional[CoreServiceRpc] = None
        session_id: Optional[str] = None
        with self._lock:
            if self._state == _CLOSED:
                error = CoreGoError(Stage.LIFECYCLE, ReasonCode.CLOSED, "core is closed; refusing to upsert")
            elif self._state != _RUNNING or self._session_id is None:
                error = CoreGoError(Stage.LIFECYCLE, ReasonCode.NOT_STARTED, "run is not running; refusing to upsert")
            else:
                error = None
                rpc = self._rpc
                session_id = self._session_id
            warn = error is not None and self._warn_once(name)
            if error is not None:
                self._last_error = error

        if error is not None:
            if warn:
                console.warning(f"Failed to upsert {name}: {error}")
            return

        assert rpc is not None
        assert session_id is not None
        try:
            call(rpc, session_id)
        except Exception as exc:
            failure = self._as_core_error(exc, Stage.RPC)
            with self._lock:
                self._last_error = failure
                warn = self._warn_once(name)
            if warn:
                console.warning(f"Failed to upsert {name}: {failure}")

    def _warn_once(self, name: str) -> bool:
        """调用方需持锁；同一方法只 warning 一次。"""
        if name in self._warned:
            return False
        self._warned.add(name)
        return True

    # ---------------------------------- 通知运行结束 ----------------------------------

    def _finish_when_local(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        return self._reject_finish("local")

    def _finish_when_offline(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        return self._reject_finish("offline")

    def _finish_when_online(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        with self._lock:
            if self._state != _RUNNING or self._session_id is None:
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
            if self._state == _RUNNING and self._session_id == session_id:
                self._state = _FINISHED
        return response

    def _reject_finish(self, mode: str) -> DeliverRunFinishResponse:
        error = CoreGoError(
            Stage.LIFECYCLE,
            ReasonCode.MODE_UNSUPPORTED,
            f"core_go only implements online mode in this build; mode={mode!r} is not handled by the Go core",
        )
        with self._lock:
            self._last_error = error
        return DeliverRunFinishResponse(success=False, message=str(error))

    # ---------------------------------- 进度查询 ----------------------------------

    def _get_operation_when_online(self) -> GetOperationStatsResponse:
        with self._lock:
            session_id = self._session_id
            rpc = self._rpc
            state = self._state
        if session_id is None or rpc is None:
            reason = ReasonCode.CLOSED if state == _CLOSED else ReasonCode.NOT_STARTED
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
            if self._state != _FINISHED or self._session_id is None:
                reason = ReasonCode.CLOSED if self._state == _CLOSED else ReasonCode.NOT_STARTED
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

    # ---------------------------------- 生命周期回收 ----------------------------------

    def close(self) -> None:
        """正常关闭：teardown → launcher 回收；不隐式 confirm，不抛出。"""
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (_CLOSED, _CLOSING, _ABORTING):
                return
            if self._state == _STARTING:
                abort_startup = True
                rpc = None
                launcher = None
            elif self._state == _NEW:
                self._state = _CLOSED
                return
            else:
                abort_startup = False
                self._state = _CLOSING
                rpc = self._rpc
                launcher = self._launcher

        if abort_startup:
            self.abort()
            return

        try:
            if rpc is not None:
                with safe.block(level="debug", message="core teardown failed during close"):
                    rpc.teardown(timeout=TEARDOWN_TIMEOUT)
            if launcher is not None:
                with safe.block(level="debug", message="core launcher close failed"):
                    launcher.close()
        finally:
            with self._lock:
                self._rpc = None
                self._launcher = None
                self._session_id = None
                self._state = _CLOSED

    def abort(self, cause: Optional[BaseException] = None) -> None:
        """失败路径回收：短 deadline 尽力 teardown，再由 launcher 有界回收；不抛出。"""
        if not self._owned_by_current_process():
            return

        with self._lock:
            if self._state in (_ABORTING, _CLOSING, _CLOSED):
                return
            self._state = _ABORTING
            if cause is not None:
                self._abort_cause = cause
            rpc = self._rpc
            launcher = self._launcher

        try:
            if rpc is not None:
                with safe.block(level="debug", message="core teardown failed during abort"):
                    rpc.teardown(timeout=TEARDOWN_TIMEOUT)
            if launcher is not None:
                with safe.block(level="debug", message="core launcher abort failed"):
                    launcher.abort()
        finally:
            with self._lock:
                self._rpc = None
                self._launcher = None
                self._session_id = None
                self._state = _CLOSED

    def _owned_by_current_process(self) -> bool:
        # 只拦真正的 os.fork()；macOS 默认 spawn，不会继承本对象。
        return not fork.is_forked(self._owner_pid)

    # ---------------------------------- 进程 fork ----------------------------------

    def fork(self) -> "CoreProtocol":
        # 不取锁：fork 时锁可能由子进程不存在的线程持有。
        return _DetachedCoreGo(cast(ModeType, self._mode), self._binary)

    @staticmethod
    def _as_core_error(exc: BaseException, stage: Stage) -> CoreGoError:
        if isinstance(exc, CoreGoError):
            return exc
        return CoreGoError(stage, ReasonCode.GRPC_ERROR, str(exc), cause=exc)
