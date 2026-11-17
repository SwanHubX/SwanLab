"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:30
@description: SwanLab 核心模块
"""

import threading
from typing import TYPE_CHECKING

from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.protocol import CoreProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext

__all__ = ["create_core", "reset_fallback_warning"]

_fallback_warn_lock = threading.Lock()
_fallback_warned = False


def reset_fallback_warning() -> None:
    """重置进程级回退告警标记（供测试隔离使用）。"""
    global _fallback_warned
    with _fallback_warn_lock:
        _fallback_warned = False


def _warn_fallback_once(reason: BaseException) -> None:
    """Go core 二进制缺失时仅告警一次，避免同进程多次 init 重复刷屏。"""
    global _fallback_warned
    with _fallback_warn_lock:
        if _fallback_warned:
            return
        _fallback_warned = True
    console.warning(f"Go core is unavailable ({reason}); falling back to the built-in Python core.")


def create_core(ctx: "RunContext") -> CoreProtocol:
    """core 对象工厂；Go core 二进制不可用时回退内置 CorePython 并告警一次。

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    mode = ctx.config.settings.mode
    if helper.get_core_impl() == "python":
        from swanlab.sdk.internal.core_python import CorePython

        return CorePython(mode=mode)

    from swanlab.sdk.internal.core_go import CoreGo
    from swanlab.sdk.internal.core_go.config import DETACH_IDLE_TIMEOUT, CoreGoError
    from swanlab.sdk.internal.core_go.runtime import locate_core_binary
    from swanlab.sdk.internal.settings.core import CoreSettings

    # Go core 是 CorePython 的全模式替代；disabled 不落盘也不起进程，无需二进制。
    if mode == "disabled":
        binary = None
    else:
        try:
            binary = locate_core_binary()
        except CoreGoError as exc:
            from swanlab.sdk.internal.core_python import CorePython

            _warn_fallback_once(exc)
            return CorePython(mode=mode)

    core_settings = getattr(ctx.config.settings, "core", None)
    detach = core_settings.detach if isinstance(core_settings, CoreSettings) else True
    idle_timeout = core_settings.idle_timeout if isinstance(core_settings, CoreSettings) else DETACH_IDLE_TIMEOUT
    return CoreGo(mode=mode, binary=binary, detach=detach, idle_timeout=idle_timeout)
