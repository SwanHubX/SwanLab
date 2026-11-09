"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:30
@description: SwanLab 核心模块
"""

from typing import TYPE_CHECKING

from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.protocol import CoreProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


def create_core(ctx: "RunContext") -> CoreProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    mode = ctx.config.settings.mode
    if helper.get_core_impl() == "python":
        from swanlab.sdk.internal.core_python import CorePython

        return CorePython(mode=mode)

    from swanlab.sdk.internal.core_go import CoreGo
    from swanlab.sdk.internal.core_go.locate import locate_core_binary

    # Go core 是 CorePython 的全模式替代；disabled 不落盘也不起进程，无需二进制。
    binary = None if mode == "disabled" else locate_core_binary()
    return CoreGo(mode=mode, binary=binary)
