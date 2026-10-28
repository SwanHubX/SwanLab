"""
@author: Nexisato
@file: core_service.py
@time: 2026/10/27 15:40
@description: CoreService 的 RPC 封装。

构造 protobuf 请求并复用 Channel 的 deadline 与错误映射。
所有 run 级方法显式携带 session_id；每个方法只调用一次 stub，
不重试、不重连，success=False 在这里转成协议错误。
"""

from typing import List

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishRequest,
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsRequest,
    GetOperationStatsResponse,
    SpinupServiceRequest,
    TeardownServiceRequest,
    UpsertColumnsRequest,
    UpsertLogsRequest,
    UpsertMediaRequest,
    UpsertSavesRequest,
    UpsertScalarsRequest,
)
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.core_go.config import TEARDOWN_TIMEOUT
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.typings.core_go import RpcChannelType

__all__ = ["CoreServiceRpc"]


class CoreServiceRpc:
    """CoreService 的强类型薄封装，不持有 channel 生命周期。"""

    def __init__(self, channel: RpcChannelType) -> None:
        # core_pb2_grpc 依赖 grpcio，延迟到真正构造 RPC 客户端时再导入。
        from swanlab.proto.swanlab.grpc.core.v1.core_pb2_grpc import CoreServiceStub

        self._channel = channel
        self._stub = CoreServiceStub(channel.raw)

    def spinup(self) -> None:
        """把服务从 NOT_READY 推进到 READY；READY 下幂等。"""
        self._channel.call(
            self._stub.SpinupService,
            SpinupServiceRequest(),
            stage=Stage.SPINUP,
        )

    def teardown(self, timeout: float = TEARDOWN_TIMEOUT) -> None:
        """请求关闭整个 core；只发请求，不等进程退出（回收由 launcher 负责）。"""
        self._channel.call(
            self._stub.TeardownService,
            TeardownServiceRequest(),
            timeout=timeout,
            stage=Stage.TEARDOWN,
        )

    def deliver_run_start(self, request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        """创建 run session，返回包含非空 session_id 的响应。"""
        response = self._channel.call(
            self._stub.DeliverRunStart,
            request,
            stage=Stage.RUN_START,
        )
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
        return response

    def upsert_columns(self, columns: List[ColumnRecord], session_id: str) -> None:
        self._channel.call(
            self._stub.UpsertColumns,
            UpsertColumnsRequest(columns=columns, session_id=session_id),
        )

    def upsert_scalars(self, scalars: List[ScalarRecord], session_id: str) -> None:
        self._channel.call(
            self._stub.UpsertScalars,
            UpsertScalarsRequest(data=scalars, session_id=session_id),
        )

    def upsert_media(self, media: List[MediaRecord], session_id: str) -> None:
        self._channel.call(
            self._stub.UpsertMedia,
            UpsertMediaRequest(data=media, session_id=session_id),
        )

    def upsert_logs(self, logs: List[LogRecord], session_id: str) -> None:
        self._channel.call(
            self._stub.UpsertLogs,
            UpsertLogsRequest(logs=logs, session_id=session_id),
        )

    def upsert_saves(self, saves: List[SaveRecord], session_id: str) -> None:
        self._channel.call(
            self._stub.UpsertSaves,
            UpsertSavesRequest(saves=saves, session_id=session_id),
        )

    def deliver_run_finish(
        self,
        request: DeliverRunFinishRequest,
        session_id: str,
    ) -> DeliverRunFinishResponse:
        """交付 run 结束事件；调用方随后用 stats 轮询到 FINISHED。"""
        request.session_id = session_id
        response = self._channel.call(
            self._stub.DeliverRunFinish,
            request,
            stage=Stage.RPC,
        )
        if not response.success:
            raise CoreGoError(
                Stage.RPC,
                ReasonCode.RUN_FINISH_REJECTED,
                response.message or "DeliverRunFinish rejected",
            )
        return response

    def get_operation_stats(self, session_id: str) -> GetOperationStatsResponse:
        response = self._channel.call(
            self._stub.GetOperationStats,
            GetOperationStatsRequest(session_id=session_id),
        )
        if not response.success:
            raise CoreGoError(
                Stage.RPC,
                ReasonCode.PROTOCOL_VIOLATION,
                response.message or "GetOperationStats rejected",
            )
        return response

    def confirm_run_finish(self, session_id: str) -> ConfirmRunFinishResponse:
        """释放 run session；不关闭 gRPC server。"""
        response = self._channel.call(
            self._stub.ConfirmRunFinish,
            ConfirmRunFinishRequest(session_id=session_id),
        )
        if not response.success:
            raise CoreGoError(
                Stage.RPC,
                ReasonCode.RUN_FINISH_REJECTED,
                response.message or "ConfirmRunFinish rejected",
            )
        return response
