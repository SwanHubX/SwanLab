"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: gRPC 边界。

Channel 是唯一 import grpc 的地方：延迟导入、channel-ready 等待、
单次调用附加 deadline，并把 grpc 状态码映射为结构化错误。
消息上限与默认超时取自 config.py，与 Go 侧保持一致。
不持有 Popen，不做任何文件系统操作。
"""

import threading
from typing import Any, Callable

from swanlab.sdk.internal.core_go.config import CHANNEL_READY_TIMEOUT, GRPC_MAX_MESSAGE_BYTES, RPC_TIMEOUT
from swanlab.sdk.internal.core_go.discovery import Endpoint
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage, reason_for_grpc

__all__ = ["Channel"]


class Channel:
    """对 grpc.Channel 的薄封装：延迟导入、ready 等待、单次调用与错误映射。"""

    def __init__(self, raw: Any, grpc_module: Any, target: str) -> None:
        self._raw = raw
        self._grpc = grpc_module
        self._target = target
        self._lock = threading.Lock()
        self._closed = False

    @classmethod
    def open(cls, endpoint: Endpoint, timeout: float = CHANNEL_READY_TIMEOUT) -> "Channel":
        """建立 channel 并等待 READY；超时或依赖缺失都转成 CoreGoError。"""
        try:
            import grpc
        except ImportError as exc:
            raise CoreGoError(
                Stage.CHANNEL,
                ReasonCode.GRPC_UNAVAILABLE,
                "grpcio is required for CoreGo online mode",
                cause=exc,
            ) from exc

        target = endpoint.target
        options = [
            ("grpc.max_receive_message_length", GRPC_MAX_MESSAGE_BYTES),
            ("grpc.max_send_message_length", GRPC_MAX_MESSAGE_BYTES),
        ]
        try:
            raw = grpc.insecure_channel(target, options=options)
        except Exception as exc:  # noqa: BLE001 - 统一转成结构化错误
            raise CoreGoError(
                Stage.CHANNEL,
                ReasonCode.GRPC_UNAVAILABLE,
                f"failed to create grpc channel for {target}: {exc}",
                cause=exc,
            ) from exc

        try:
            grpc.channel_ready_future(raw).result(timeout=timeout)
        except grpc.FutureTimeoutError as exc:
            raw.close()
            raise CoreGoError(
                Stage.CHANNEL,
                ReasonCode.CHANNEL_TIMEOUT,
                f"grpc channel not ready within {timeout:.1f}s: {target}",
                cause=exc,
            ) from exc
        except Exception as exc:  # noqa: BLE001 - ready 阶段的错误统一归类
            raw.close()
            raise CoreGoError(
                Stage.CHANNEL,
                ReasonCode.GRPC_UNAVAILABLE,
                f"grpc channel failed before ready: {target}: {exc}",
                cause=exc,
            ) from exc

        return cls(raw, grpc, target)

    @property
    def raw(self) -> Any:
        return self._raw

    def call(
        self,
        method: Callable[..., Any],
        request: Any,
        timeout: float = RPC_TIMEOUT,
        *,
        stage: Stage = Stage.RPC,
    ) -> Any:
        """调用一次 stub 方法；不重试、不重连，gRPC 错误只映射一次。"""
        if self._closed:
            raise CoreGoError(stage, ReasonCode.CLOSED, f"grpc channel is closed: {self._target}")

        try:
            return method(request, timeout=timeout)
        except self._grpc.RpcError as exc:
            raise self._map_rpc_error(exc, stage) from exc

    def _map_rpc_error(self, exc: Any, stage: Stage) -> CoreGoError:
        code = exc.code()
        name = getattr(code, "name", None) or "UNKNOWN"
        value = getattr(code, "value", None)
        grpc_code = value[0] if isinstance(value, tuple) and value else None
        details = exc.details() if hasattr(exc, "details") else None
        return CoreGoError(
            stage,
            reason_for_grpc(name),
            details or f"gRPC call failed with {name}",
            grpc_name=name,
            grpc_code=grpc_code,
        )

    def close(self) -> None:
        """幂等关闭；关闭失败不抛出。"""
        with self._lock:
            if self._closed:
                return
            self._closed = True
            raw = self._raw
        try:
            raw.close()
        except Exception:  # noqa: BLE001 - 关闭失败不影响回收流程
            pass
