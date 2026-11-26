"""
@author: Nexisato
@file: _common.py
@time: 2026/11/26
@description: core_go RPC 服务封装共用的带状态校验调用辅助。
"""

from typing import Any, Callable, Optional, Protocol, TypeVar

from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.typings.core_go import RpcChannelType


class _StatusResponse(Protocol):
    @property
    def success(self) -> bool: ...

    @property
    def message(self) -> str: ...


_RespT = TypeVar("_RespT", bound=_StatusResponse)


def call_checked(
    channel: RpcChannelType,
    method: Callable[..., _RespT],
    request: Any,
    *,
    timeout: Optional[float],
    stage: Stage = Stage.RPC,
    reject_reason: ReasonCode,
    default_message: str,
) -> _RespT:
    """执行一次 RPC 调用并校验 response.success；失败时统一转为 CoreGoError。"""
    response = channel.call(method, request, timeout=timeout, stage=stage)
    if not response.success:
        raise CoreGoError(stage, reject_reason, response.message or default_message)
    return response
