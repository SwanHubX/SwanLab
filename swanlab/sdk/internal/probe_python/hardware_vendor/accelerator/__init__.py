"""
@author: cunyue
@file: __init__.py
@time: 2026/3/31 01:49
@description: 计算加速卡信息模块
"""

from typing import Dict, Type

from swanlab.sdk.typings.probe_python import AcceleratorVendor
from swanlab.sdk.typings.probe_python.hardware_vendor import AcceleratorProtocol

from .amd import AMDGPU
from .cambricon import CambriconMLU
from .huawei import AscendNPU
from .hygon import HygonDCU
from .iluvatar import IluvatarGPU
from .kunlunxin import KunlunxinXPU
from .metax import MetaXGPU
from .moorethreads import MooreThreadsGPU
from .nvidia import NvidiaGPU

__all__ = ["ACCELERATOR_REGISTRY"]

ACCELERATOR_REGISTRY: Dict[AcceleratorVendor, Type[AcceleratorProtocol]] = {
    "nvidia": NvidiaGPU,
    "rocm": AMDGPU,
    "ascend": AscendNPU,
    "cambricon": CambriconMLU,
    "iluvatar": IluvatarGPU,
    "metax": MetaXGPU,
    "moorethreads": MooreThreadsGPU,
    "kunlunxin": KunlunxinXPU,
    "hygon": HygonDCU,
}
