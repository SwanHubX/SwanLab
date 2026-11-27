"""
@author: Nexisato
@file: channel.py
@time: 2026/10/27 15:40
@description: gRPC 边界。

Channel 是唯一 import grpc 的地方：延迟导入、channel-ready 等待、
单次调用附加 deadline，并把 grpc 状态码映射为结构化错误。
消息上限与默认超时取自 config，与 Go 侧保持一致。
不持有 Popen，不做任何文件系统操作。
"""

import os
import threading
from typing import Any, Callable, Optional

from swanlab.sdk.internal.core_go.config import (
    CHANNEL_READY_TIMEOUT,
    GRPC_MAX_MESSAGE_BYTES,
    RPC_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
    reason_for_grpc,
)
from swanlab.sdk.internal.core_go.runtime import Endpoint
from swanlab.sdk.internal.pkg import safe

__all__ = ["Channel"]

# grpc C-core 默认开启 fork 支持，fork 时可能向 stderr 输出 absl 日志；
# SwanLab 的 fork 子进程禁止 Core RPC，这里在首次导入 grpc 前默认关闭，用户显式设置优先。
os.environ.setdefault("GRPC_ENABLE_FORK_SUPPORT", "0")

# gRPC 默认在空闲 30 分钟后断开连接；设置为 int32 上限（约 24.9 天）以覆盖长时训练空闲窗口。
_IDLE_TIMEOUT_DISABLED_MS = 2**31 - 1


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
                ReasonCode.GRPC_LOCAL_ERROR,
                "grpcio is required for CoreGo",
                cause=exc,
            ) from exc

        target = endpoint.target
        options = [
            ("grpc.max_receive_message_length", GRPC_MAX_MESSAGE_BYTES),
            ("grpc.max_send_message_length", GRPC_MAX_MESSAGE_BYTES),
            ("grpc.client_idle_timeout_ms", _IDLE_TIMEOUT_DISABLED_MS),
        ]
        try:
            raw = grpc.insecure_channel(target, options=options)
        except Exception as exc:  # noqa: BLE001 - 统一转成结构化错误
            # 建不出 channel 属于本机环境/参数问题，重试与换进程都无意义
            raise CoreGoError(
                Stage.CHANNEL,
                ReasonCode.GRPC_LOCAL_ERROR,
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
        timeout: Optional[float] = RPC_TIMEOUT,
        *,
        stage: Stage = Stage.RPC,
    ) -> Any:
        """调用一次 stub 方法；不重试、不重连，gRPC 错误只映射一次。"""
        if self._closed:
            raise CoreGoError(stage, ReasonCode.CLOSED, f"grpc channel is closed: {self._target}")

        self._check_request_size(request, stage)
        try:
            return method(request, timeout=timeout)
        except self._grpc.RpcError as exc:
            raise self._map_rpc_error(exc, stage) from exc

    @staticmethod
    def _check_request_size(request: Any, stage: Stage) -> None:
        """发送前拦截超过 gRPC 上限的请求，避免进入重试队列。"""
        byte_size = getattr(request, "ByteSize", None)
        if byte_size is None:
            return
        size = byte_size()
        if size > GRPC_MAX_MESSAGE_BYTES:
            raise CoreGoError(
                stage,
                ReasonCode.REQUEST_TOO_LARGE,
                f"request size {size} exceeds gRPC limit {GRPC_MAX_MESSAGE_BYTES}",
            )

    def _map_rpc_error(self, exc: Any, stage: Stage) -> CoreGoError:
        code = exc.code()
        name = getattr(code, "name", None) or "UNKNOWN"
        value = getattr(code, "value", None)
        grpc_code = value[0] if isinstance(value, tuple) and value else None
        details = exc.details() if hasattr(exc, "details") else None
        reason = reason_for_grpc(name)
        # 客户端消息超限也报 RESOURCE_EXHAUSTED，单独归为超限：上层可以拆片重发
        if reason is ReasonCode.RESOURCE_EXHAUSTED and "larger than max" in (details or "").lower():
            reason = ReasonCode.REQUEST_TOO_LARGE
        return CoreGoError(
            stage,
            reason,
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
        with safe.block(level="debug", message=f"failed to close grpc channel {self._target}"):
            raw.close()
