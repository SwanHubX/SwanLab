"""
@author: cunyue
@file: hardware_vendor.py
@time: 2026/3/31 13:14
@description: SwanLab 运行时系统组件 - 硬件供应商相关的类型定义和接口规范
"""

from abc import ABC, abstractmethod
from typing import Optional

from swanlab.sdk.typings.run.system import AcceleratorSnapshot, AppleSiliconSnapshot, CPUSnapshot, MemorySnapshot


class CpuProtocol(ABC):
    """
    CPU 相关的接口规范
    """

    @abstractmethod
    def get(self) -> Optional[CPUSnapshot]: ...


class MemoryProtocol(ABC):
    """
    内存相关的接口规范
    """

    @abstractmethod
    def get(self) -> Optional[MemorySnapshot]: ...


class AppleSiliconProtocol(ABC):
    """
    苹果统一芯片（Apple Silicon）相关的接口规范
    """

    @abstractmethod
    def get(self) -> Optional[AppleSiliconSnapshot]: ...


class AcceleratorProtocol(ABC):
    """
    加速器相关的接口规范
    """

    @abstractmethod
    def get(self) -> Optional[AcceleratorSnapshot]: ...
