"""
@author: cunyue
@file: __init__.py
@time: 2026/4/21 15:25
@description: SwanLab 系统监控模块
"""

from typing import Optional

from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.protocol import CoreProtocol, ProbeProtocol


def create_probe(core: Optional[CoreProtocol] = None, disabled: bool = False) -> ProbeProtocol:
    """core对象工厂

    :param ctx: 运行上下文，包含配置信息和运行时状态
    """
    if helper.get_probe_impl() == "python":
        from swanlab.sdk.internal.probe_python import ProbePython

        return ProbePython(core, disabled)
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError("The SwanLab Rust probe runtime is not available yet.")
