"""
@author: Nexisato
@file: sync_service.py
@time: 2026/11/18
@description: CoreSyncService 的 RPC 封装。

构造 protobuf 请求并复用 Channel 的 deadline 与错误映射。
服务端按 gRPC 连接绑定唯一 sync 会话，Flush/Stats/Confirm 均携带 google.protobuf.Empty；
每个方法只调用一次 stub，不重试、不重连，success=False 在这里转成协议错误。
"""

from google.protobuf.empty_pb2 import Empty

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
    DeliverSyncStartResponse,
)
from swanlab.sdk.internal.core_go.config import (
    CONFIRM_TIMEOUT,
    RUN_START_TIMEOUT,
    STATS_TIMEOUT,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.rpc._common import call_checked
from swanlab.sdk.typings.core_go import RpcChannelType

__all__ = ["CoreSyncRpc"]


class CoreSyncRpc:
    """CoreSyncService 的强类型薄封装，不持有 channel 生命周期。"""

    def __init__(self, channel: RpcChannelType) -> None:
        # sync_pb2_grpc 依赖 grpcio，延迟到真正构造 RPC 客户端时再导入。
        from swanlab.proto.swanlab.grpc.core.v1.sync_pb2_grpc import CoreSyncServiceStub

        self._channel = channel
        self._stub = CoreSyncServiceStub(channel.raw)

    def deliver_sync_start(self, request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        """校验本地配置与 run-*.swanlab 文件头并创建连接绑定的 sync 会话。"""
        return call_checked(
            self._channel,
            self._stub.DeliverSyncStart,
            request,
            timeout=RUN_START_TIMEOUT,
            stage=Stage.RUN_START,
            reject_reason=ReasonCode.RUN_START_REJECTED,
            default_message="DeliverSyncStart rejected",
        )

    def deliver_sync_flush(self) -> DeliverSyncFlushResponse:
        """完成云端实验创建/恢复与摘要水位初始化，并启动后台扫描上传。"""
        return call_checked(
            self._channel,
            self._stub.DeliverSyncFlush,
            Empty(),
            timeout=RUN_START_TIMEOUT,
            stage=Stage.RUN_START,
            reject_reason=ReasonCode.RUN_START_REJECTED,
            default_message="DeliverSyncFlush rejected",
        )

    def get_operation_stats(self) -> GetOperationStatsResponse:
        """查询当前连接绑定 sync 会话的上传进度与运行状态快照。"""
        return call_checked(
            self._channel,
            self._stub.GetOperationStats,
            Empty(),
            timeout=STATS_TIMEOUT,
            reject_reason=ReasonCode.PROTOCOL_VIOLATION,
            default_message="GetOperationStats rejected",
        )

    def confirm_sync_finish(self) -> ConfirmSyncFinishResponse:
        """等待扫描与上传全部排空并上报终态（遵循 CONFIRM_TIMEOUT=None 无限等待语义）。"""
        return call_checked(
            self._channel,
            self._stub.ConfirmSyncFinish,
            Empty(),
            timeout=CONFIRM_TIMEOUT,
            reject_reason=ReasonCode.RUN_FINISH_REJECTED,
            default_message="ConfirmSyncFinish rejected",
        )
