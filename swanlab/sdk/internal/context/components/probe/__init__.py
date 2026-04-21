"""
@author: cunyue
@file: __init__.py
@time: 2026/4/21 15:25
@description: SwanLab 系统监控模块
"""

from typing import TYPE_CHECKING

from swanlab.sdk.protocol import ProbeEnum, ProbeProtocol

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


# TODO: 未来实现probe以后，python版本依旧会有一段时间的同时存在时间。后续实现一种机制，选择不同的probe实现
probe_enum: ProbeEnum = ProbeEnum.PROBE_PYTHON


def create_probe(ctx: "RunContext") -> ProbeProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if probe_enum == ProbeEnum.PROBE_PYTHON:
        from swanlab.sdk.internal.probe_python import ProbePython

        return ProbePython(ctx)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError(f"CoreEnum {probe_enum} is not supported yet.")
