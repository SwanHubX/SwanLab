"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27
@description: core_go RPC 边界的可复用类型提示。

只描述 CoreServiceRpc 实际依赖的 channel 能力，不绑定具体 grpc 实现，
真实 Channel 与测试替身都可以满足该协议。
"""

from typing import Any, Callable, Protocol

__all__ = ["RpcChannelType"]


class RpcChannelType(Protocol):
    """CoreServiceRpc 需要的最小 channel 能力。"""

    @property
    def raw(self) -> Any: ...

    def call(
        self,
        method: Callable[..., Any],
        request: Any,
        timeout: float = ...,
        *,
        # stage 是 core_go 内部的 Stage 枚举；typings 不反向依赖 internal，
        # 这里保留 Any，真实 Channel 会把它收窄为 Stage。
        stage: Any = ...,
    ) -> Any: ...
