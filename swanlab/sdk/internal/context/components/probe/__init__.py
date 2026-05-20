"""
@author: cunyue
@file: __init__.py
@time: 2026/4/21 15:25
@description: SwanLab 系统监控模块
"""

from typing import TYPE_CHECKING

from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.protocol import ProbeProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


def create_probe(ctx: "RunContext") -> ProbeProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if helper.get_probe_impl() == "python":
        from swanlab.sdk.internal.probe_python import ProbePython

        return ProbePython(ctx.config.settings.mode, ctx.core)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError("The SwanLab Rust probe runtime is not available yet.")
