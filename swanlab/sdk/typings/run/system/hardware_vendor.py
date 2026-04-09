"""
@author: cunyue
@file: hardware_vendor.py
@time: 2026/3/31 13:14
@description: SwanLab 运行时系统组件 - 硬件供应商相关的类型定义和接口规范
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

from swanlab.sdk.typings.run.system import (
    AcceleratorSnapshot,
    AppleSiliconSnapshot,
    CPUSnapshot,
    MemorySnapshot,
    SystemShim,
)

if TYPE_CHECKING:
    from swanlab.sdk import Run


class CpuProtocol(ABC):
    """
    CPU 相关的接口规范
    """

    def __init__(self, shim: SystemShim):
        self._shim = shim

    @classmethod
    @abstractmethod
    def new(cls, run: "Run", shim: SystemShim) -> Optional["CpuProtocol"]: ...

    @abstractmethod
    def collect(self): ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[CPUSnapshot]: ...


class MemoryProtocol(ABC):
    """
    内存相关的接口规范
    """

    def __init__(self, shim: SystemShim):
        self._shim = shim

    @classmethod
    @abstractmethod
    def new(cls, run: "Run", shim: SystemShim) -> Optional["MemoryProtocol"]: ...

    @abstractmethod
    def collect(self): ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[MemorySnapshot]: ...


class AppleSiliconProtocol(ABC):
    """
    苹果统一芯片（Apple Silicon）相关的接口规范
    """

    def __init__(self, shim: SystemShim):
        self._shim = shim

    @classmethod
    @abstractmethod
    def new(cls, run: "Run", shim: SystemShim) -> Optional["AppleSiliconProtocol"]: ...

    @abstractmethod
    def collect(self): ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[AppleSiliconSnapshot]: ...


class AcceleratorProtocol(ABC):
    """
    加速器相关的接口规范
    """

    def __init__(self, shim: SystemShim):
        self._shim = shim

    @classmethod
    @abstractmethod
    def new(cls, run: "Run", shim: SystemShim) -> Optional["AcceleratorProtocol"]: ...

    @abstractmethod
    def collect(self): ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[AcceleratorSnapshot]: ...
