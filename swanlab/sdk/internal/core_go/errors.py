"""
@author: Nexisato
@file: errors.py
@time: 2026/10/27 15:40
@description: 结构化错误与原因码。

CoreGoError 携带阶段（stage）、原因码（reason）与可选的 gRPC 状态码，
让上层能区分启动失败、通道超时、RPC 拒绝与 UNIMPLEMENTED。
内部专用，不复用 swanlab/exceptions.py 中的 HTTP/DataStore 异常。
"""

from enum import Enum
from typing import Any, Dict, Optional

__all__ = ["CoreGoError", "ReasonCode", "Stage", "reason_for_grpc"]


class Stage(str, Enum):
    """失败发生的阶段，用于快速定位握手链路中的位置。"""

    RUNTIME = "RUNTIME"
    SPAWN = "SPAWN"
    PORT_FILE = "PORT_FILE"
    CHANNEL = "CHANNEL"
    SPINUP = "SPINUP"
    RUN_START = "RUN_START"
    RPC = "RPC"
    TEARDOWN = "TEARDOWN"
    REAP = "REAP"
    CLEANUP = "CLEANUP"
    FORK = "FORK"
    LIFECYCLE = "LIFECYCLE"


class ReasonCode(str, Enum):
    """机器可判别的原因码，同一个原因在不同阶段含义一致。"""

    # 启动
    BINARY_NOT_FOUND = "BINARY_NOT_FOUND"
    RUNTIME_UNAVAILABLE = "RUNTIME_UNAVAILABLE"
    SPAWN_FAILED = "SPAWN_FAILED"
    EARLY_EXIT = "EARLY_EXIT"
    ENDPOINT_BUSY = "ENDPOINT_BUSY"
    USAGE_ERROR = "USAGE_ERROR"
    # 发现
    PORT_FILE_TIMEOUT = "PORT_FILE_TIMEOUT"
    PORT_FILE_INVALID = "PORT_FILE_INVALID"
    OWNERSHIP_MISMATCH = "OWNERSHIP_MISMATCH"
    # 连接
    GRPC_UNAVAILABLE = "GRPC_UNAVAILABLE"
    CHANNEL_TIMEOUT = "CHANNEL_TIMEOUT"
    # RPC
    DEADLINE_EXCEEDED = "DEADLINE_EXCEEDED"
    UNIMPLEMENTED = "UNIMPLEMENTED"
    FAILED_PRECONDITION = "FAILED_PRECONDITION"
    SESSION_NOT_FOUND = "SESSION_NOT_FOUND"
    RESOURCE_EXHAUSTED = "RESOURCE_EXHAUSTED"
    INVALID_ARGUMENT = "INVALID_ARGUMENT"
    GRPC_ERROR = "GRPC_ERROR"
    # 语义
    PROTOCOL_VIOLATION = "PROTOCOL_VIOLATION"
    RUN_START_REJECTED = "RUN_START_REJECTED"
    RUN_FINISH_REJECTED = "RUN_FINISH_REJECTED"
    # 客户端状态
    NOT_STARTED = "NOT_STARTED"
    ALREADY_STARTED = "ALREADY_STARTED"
    CLOSED = "CLOSED"
    MODE_UNSUPPORTED = "MODE_UNSUPPORTED"
    FORK_UNSUPPORTED = "FORK_UNSUPPORTED"


# gRPC 状态名到原因码的唯一映射，模块内不 import grpc。
_GRPC_REASON = {
    "UNIMPLEMENTED": ReasonCode.UNIMPLEMENTED,
    "UNAVAILABLE": ReasonCode.GRPC_UNAVAILABLE,
    "DEADLINE_EXCEEDED": ReasonCode.DEADLINE_EXCEEDED,
    "FAILED_PRECONDITION": ReasonCode.FAILED_PRECONDITION,
    "NOT_FOUND": ReasonCode.SESSION_NOT_FOUND,
    "RESOURCE_EXHAUSTED": ReasonCode.RESOURCE_EXHAUSTED,
    "INVALID_ARGUMENT": ReasonCode.INVALID_ARGUMENT,
}


def reason_for_grpc(name: str) -> ReasonCode:
    """把 gRPC 状态名转换为 ReasonCode；未知状态统一落到 GRPC_ERROR。

    纯字符串映射，单元测试不需要安装 grpcio。
    """
    return _GRPC_REASON.get(str(name).strip().upper(), ReasonCode.GRPC_ERROR)


class CoreGoError(Exception):
    """CoreGo 握手链路的唯一异常类型。

    Attributes:
        stage: 失败阶段。
        reason: 机器可判别的原因码。
        message: 给人看的说明。
        grpc_name: gRPC 状态名，例如 ``"UNIMPLEMENTED"``。
        grpc_code: gRPC 规范数字码，便于跨语言对照。
        detail: 有界诊断信息，例如 stderr 尾部。
    """

    def __init__(
        self,
        stage: Stage,
        reason: ReasonCode,
        message: str,
        *,
        grpc_name: Optional[str] = None,
        grpc_code: Optional[int] = None,
        detail: Optional[str] = None,
        cause: Optional[BaseException] = None,
    ) -> None:
        super().__init__(message)
        self.stage = stage
        self.reason = reason
        self.message = message
        self.grpc_name = grpc_name
        self.grpc_code = grpc_code
        self.detail = detail
        if cause is not None:
            self.__cause__ = cause

    def __str__(self) -> str:
        text = f"[{self.stage.value}/{self.reason.value}] {self.message}"
        if self.grpc_name:
            text += f" (grpc={self.grpc_name})"
        return text

    def to_dict(self) -> Dict[str, Any]:
        """返回可 JSON 序列化的快照，供 harness 与日志使用。"""
        return {
            "stage": self.stage.value,
            "reason": self.reason.value,
            "message": self.message,
            "grpc_name": self.grpc_name,
            "grpc_code": self.grpc_code,
            "detail": self.detail,
        }
