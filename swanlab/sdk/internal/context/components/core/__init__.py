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
    if helper.get_core_impl() == "python":
        from swanlab.sdk.internal.core_python import CorePython

        return CorePython(mode=ctx.config.settings.mode)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError("The SwanLab Go core runtime is not available yet.")
