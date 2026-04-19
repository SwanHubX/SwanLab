"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 22:30
@description: SwanLab 核心模块
"""

from typing import TYPE_CHECKING

from swanlab.sdk.protocol import CoreEnum, CoreProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


# TODO: 未来实现core以后，python版本依旧会有一段时间的同时存在时间。后续实现一种机制，选择不同的core实现
core_enum: CoreEnum = CoreEnum.CORE_PYTHON


def create_core(ctx: "RunContext") -> CoreProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if core_enum == CoreEnum.CORE_PYTHON:
        from swanlab.sdk.internal.core_python import CorePython

        return CorePython(ctx)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError(f"CoreEnum {core_enum} is not supported yet.")
