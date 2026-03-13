"""
@author: cunyue
@file: core.py
@time: 2026/3/13
@description: SwanLab Core 接口协议，负责对产出的Record做持久化处理和后端交互
"""

from abc import ABC, abstractmethod
from typing import List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunContext

__all__ = ["CoreProtocol"]


class CoreProtocol(ABC):
    """
    SwanLab Core 的抽象基类。

    实现方负责：
    - 将 Record 批次持久化到本地（loglevel 格式）
    - 触发回调、向云端上传等后续业务
    - 生命周期管理（startup / shutdown 配对，防止重复初始化）

    调用方（BackgroundConsumer）只依赖此抽象类，不依赖具体实现。
    """

    def __init__(self, ctx: RunContext):
        self._ctx = ctx

    @abstractmethod
    def startup(
        self,
        cloud: bool,
        persistence: bool,
    ) -> None:
        """初始化资源，与 shutdown() 配对。重复调用应抛出错误。"""
        ...

    @abstractmethod
    def handle_records(self, records: List[Record]) -> None:
        """处理一批 Record：持久化、回调、上传等，调用方保证列表非空。"""
        ...

    @abstractmethod
    def shutdown(self) -> None:
        """刷盘、触发收尾回调、释放所有资源。调用后此实例不应再被使用。"""
        ...
