"""
@author: cunyue
@file: impl.py
@time: 2026/5/17
@description: SwanLab SDK component implementation helpers
"""

from swanlab.sdk.protocol import CoreEnum, ProbeEnum

# TODO: 未来实现 probe 后，Python 版本依旧会有一段时间的共存期。后续实现一种机制，选择不同的 probe 实现
_probe_impl: ProbeEnum = ProbeEnum.PROBE_PYTHON

# TODO: 未来实现 core 后，Python 版本依旧会有一段时间的共存期。后续实现一种机制，选择不同的 core 实现
_core_impl: CoreEnum = CoreEnum.CORE_PYTHON


def get_probe_impl() -> ProbeEnum:
    """获取当前 probe 实现类型。"""
    return _probe_impl


def get_core_impl() -> CoreEnum:
    """获取当前 core 实现类型。"""
    return _core_impl
