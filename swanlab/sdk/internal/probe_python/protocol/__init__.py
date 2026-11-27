"""
@author: cunyue
@file: hardware_vendor.py
@time: 2026/3/31 13:14
@description: SwanLab 运行时系统组件 - 硬件供应商相关的类型定义和接口规范
"""

from abc import ABC, abstractmethod
from typing import Callable, List, Optional, Set, Tuple, Union

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.probe_python.typings import (
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
        # 已经告警过的指标 key：同一个指标后续失败只写 debug，避免刷屏
        self._warned_failures: Set[str] = set()

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        newly_failed: List[Tuple[str, BaseException]] = []
        for key, handler in self._handlers:
            try:
                results.append((key, handler()))
            except Exception as exc:
                if key in self._warned_failures:
                    console.debug(
                        f"Failed to collect {key} while calling {handler.__name__}: {exc}", write_to_tty=False
                    )
                else:
                    self._warned_failures.add(key)
                    newly_failed.append((key, exc))
        if newly_failed:
            keys = ", ".join(key for key, _ in newly_failed)
            console.warning(f"Failed to collect system metrics: {keys} ({newly_failed[0][1]}).")
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
