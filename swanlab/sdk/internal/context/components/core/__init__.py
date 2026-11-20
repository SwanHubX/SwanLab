"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:30
@description: SwanLab 核心模块
"""

import threading
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Optional

from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.protocol import CoreProtocol
from swanlab.sdk.typings.run import ModeType

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext
    from swanlab.sdk.internal.settings import Settings

__all__ = [
    "CoreBackend",
    "create_core",
    "resolve_core_backend",
    "resolve_detach",
    "warn_fallback_once",
]

_fallback_warn_lock = threading.Lock()
_fallback_warned = False


@dataclass(frozen=True)
class CoreBackend:
    """本次运行实际使用的 core 后端；``fallback_reason`` 非空表示 Go core 不可用后回退。"""

    impl: Literal["python", "go"]
    binary: Optional[Path] = None
    fallback_reason: Optional[BaseException] = None


def resolve_core_backend(mode: ModeType) -> CoreBackend:
    """解析实际使用的 core 后端；Go core 不可用时回退内置实现。"""
    if helper.get_core_impl() == "python":
        return CoreBackend("python")
    # disabled 不落盘也不起进程，无需定位二进制。
    if mode == "disabled":
        return CoreBackend("go")

    from swanlab.sdk.internal.core_go.config import CoreGoError
    from swanlab.sdk.internal.core_go.runtime import locate_core_binary

    try:
        return CoreBackend("go", locate_core_binary())
    except CoreGoError as exc:
        return CoreBackend("python", fallback_reason=exc)


def reset_fallback_warning() -> None:
    """重置进程级回退告警标记（供测试隔离使用）。"""
    global _fallback_warned
    with _fallback_warn_lock:
        _fallback_warned = False


def resolve_detach(settings: "Settings") -> bool:
    """解析 Go core 是否需要以共享 detached daemon 运行。

    ``parallel="shared"`` 默认走 Detached Mode（多 rank 共享一个 core 进程与连接池）；
    显式设置的 ``core.detach`` 优先，便于并行场景仍选用独占 core。

    NOTE: 仅对 Go core 生效；CorePython 没有独立进程，两个字段都被静默忽略。
    """
    if settings.core.detach is not None:
        return settings.core.detach
    return settings.run.parallel == "shared"


def warn_fallback_once(reason: BaseException) -> None:
    """Go core 二进制缺失时仅告警一次，避免同进程多次 init / sync 重复刷屏。"""
    global _fallback_warned
    with _fallback_warn_lock:
        if _fallback_warned:
            return
        _fallback_warned = True
    console.warning(f"Go core is unavailable ({reason}); falling back to the built-in Python core.")


def create_core(ctx: "RunContext", backend: Optional[CoreBackend] = None) -> CoreProtocol:
    """core 对象工厂；Go core 不可用时回退内置 CorePython 并告警一次。

    :param ctx: 运行上下文，包含配置信息和运行时状态
    :param backend: 已解析的后端；未提供时按当前配置重新解析
    """
    mode = ctx.config.settings.mode
    resolved = backend if backend is not None else resolve_core_backend(mode)
    if resolved.impl == "python":
        # CorePython 跑在当前进程内，没有独立的 core 进程可共享，
        # 因此 core.detach 与 parallel="shared" 的共享语义在这里静默忽略（也不告警）。
        from swanlab.sdk.internal.core_python import CorePython

        if resolved.fallback_reason is not None:
            warn_fallback_once(resolved.fallback_reason)
        return CorePython(mode=mode)

    from swanlab.sdk.internal.core_go import CoreGo

    return CoreGo(
        mode=mode,
        binary=resolved.binary,
        detach=resolve_detach(ctx.config.settings),
        idle_timeout=ctx.config.settings.core.idle_timeout,
    )
