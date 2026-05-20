"""
@author: cunyue
@file: impl.py
@time: 2026/5/17
@description: SwanLab SDK component implementation helpers
"""

from typing import Literal

from swanlab.sdk.protocol import CoreEnum, ProbeEnum

# TODO: 未来实现 probe 后，Python 版本依旧会有一段时间的共存期。后续实现一种机制，选择不同的 probe 实现
_probe_impl: ProbeEnum = ProbeEnum.PROBE_PYTHON

# TODO: 未来实现 core 后，Python 版本依旧会有一段时间的共存期。后续实现一种机制，选择不同的 core 实现
_core_impl: CoreEnum = CoreEnum.CORE_PYTHON


def get_probe_impl() -> Literal["python", "rust"]:
    """获取当前 probe 实现类型。"""
    if _probe_impl == ProbeEnum.PROBE_PYTHON:
        return "python"
    elif _probe_impl == ProbeEnum.PROBE:
        return "rust"
    else:
        raise ValueError(f"Unknown probe impl: {_probe_impl}")


def get_core_impl() -> Literal["python", "go"]:
    """获取当前 core 实现类型。"""
    if _core_impl == CoreEnum.CORE_PYTHON:
        return "python"
    elif _core_impl == CoreEnum.CORE:
        return "go"
    else:
        raise ValueError(f"Unknown core impl: {_core_impl}")
