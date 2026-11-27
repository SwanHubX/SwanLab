"""
@author: cunyue
@file: __init__.py
@time: 2026/4/21 15:25
@description: SwanLab 系统监控模块
"""

from typing import Optional

from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.protocol import CoreProtocol, ProbeProtocol


def create_probe(core: Optional[CoreProtocol] = None, disabled: bool = False) -> ProbeProtocol:
    """probe 对象工厂；当前固定返回内置 ProbePython 实现。"""
    if helper.get_probe_impl() != "python":
        console.warning("Non-Python probe runtime is unavailable; using built-in ProbePython instead.")
    from swanlab.sdk.internal.probe_python import ProbePython

    return ProbePython(core, disabled)
