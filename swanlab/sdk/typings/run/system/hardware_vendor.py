"""
@author: cunyue
@file: hardware_vendor.py
@time: 2026/3/31 13:14
@description: SwanLab 运行时系统组件 - 硬件供应商相关的类型定义和接口规范
"""

from abc import ABC, abstractmethod
from typing import Callable, List, Optional, Tuple, Union

from swanlab.sdk.internal.pkg.safe import safe_block
from swanlab.sdk.typings.run.system import (
    AcceleratorSnapshot,
    AppleSiliconSnapshot,
    CPUSnapshot,
    MemorySnapshot,
    SystemScalars,
    SystemShim,
)

CollectResult = Tuple[str, Union[int, float]]
"""
采集结果： key, value
key 为这个指标对应的唯一标识符
"""


class CollectorProtocol(ABC):
    """
    采集器协议
    """

    def __init__(self, shim: SystemShim):
        self._shim = shim
        self._handlers: List[Tuple[str, Callable[[], Union[int, float]]]] = []

    def collect(self) -> List[CollectResult]:
        results = []
        for key, handler in self._handlers:
            with safe_block(message=f"Failed to collect while calling {handler.__name__}"):
                results.append((key, handler()))
        return results


class CpuProtocol(CollectorProtocol):
    """
    CPU 相关的接口规范
    """

    @classmethod
    @abstractmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["CpuProtocol", SystemScalars]]: ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[CPUSnapshot]: ...


class MemoryProtocol(CollectorProtocol):
    """
    内存相关的接口规范
    """

    @classmethod
    @abstractmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["MemoryProtocol", SystemScalars]]: ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[MemorySnapshot]: ...


class AppleSiliconProtocol(CollectorProtocol):
    """
    苹果统一芯片（Apple Silicon）相关的接口规范
    """

    @classmethod
    @abstractmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["AppleSiliconProtocol", SystemScalars]]: ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[AppleSiliconSnapshot]: ...


class AcceleratorProtocol(CollectorProtocol):
    """
    加速器相关的接口规范
    """

    @classmethod
    @abstractmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["AcceleratorProtocol", SystemScalars]]: ...

    @staticmethod
    @abstractmethod
    def get() -> Optional[AcceleratorSnapshot]: ...
