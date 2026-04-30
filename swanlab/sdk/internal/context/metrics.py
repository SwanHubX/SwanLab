"""
@author: cunyue
@file: metrics.py
@time: 2026/3/12 01:29
@description: SwanLab 运行时指标管理
"""

import math
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Set, Union

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.pkg import console

# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


@dataclass(**DATACLASS_KWARGS)
class BaseMetric(ABC):
    _column: ColumnRecord
    steps: Set[int] = field(default_factory=set, init=False)

    @property
    def type(self) -> ColumnType:
        return self._column.column_type

    def ensure_type_match(self, metric_type: ColumnType):
        if self.type != metric_type:
            raise TypeError(
                f"Metric '{self._column.column_key}' has already been defined as "
                f"{ColumnType.Name(self.type)}, not {ColumnType.Name(metric_type)}."
            )

    def check_and_mark_logged(self, step: int) -> bool:
        """
        确保step已经被记录，如果已经记录则返回True，否则返回False
        引入此方法是实现 https://github.com/SwanHubX/SwanLab/issues/1576 的前置条件
        我们需要等待后端准备好
        :param step: 步数
        :return: 是否已经记录
        """
        if step in self.steps:
            return True
        self.steps.add(step)
        return False

    @abstractmethod
    def update(self, data_record: Any): ...


# 使用解包的方式传入参数
@dataclass(**DATACLASS_KWARGS)
class ScalarMetric(BaseMetric):
    # 指标最新值
    latest: Optional[Union[float, int]] = None
    # 指标的最大值
    max: Optional[Union[float, int]] = None
    # 指标的最小值
    min: Optional[Union[float, int]] = None

    def update(self, data_record: ScalarRecord):
        """
        更新标量指标值
        :param data_record: 标量指标记录
        """
        value = data_record.value.number
        # 1. 如果为有效值，则更新
        if math.isfinite(value):
            self.latest = value
            if self.max is None or value > self.max:
                self.max = value
            if self.min is None or value < self.min:
                self.min = value
            return
        # 2. 如果为无效值，则忽略
        key = self._column.column_key
        console.debug(f"Invalid scalar value: {value} for metric '{key}', ignored when updating.")


@dataclass(**DATACLASS_KWARGS)
class MediaMetric(BaseMetric):
    # 媒体存储路径，绝对路径
    path: Path

    def update(self, data_record: Any):
        """
        更新媒体指标值，媒体指标目前没什么可更新的，因此直接pass
        """
        pass


# 指标状态，实验运行过程中不断更新
# 注意：此类不加锁保护。设计上对 RunMetrics 的写操作发生在以下线程：
#   - next_step: 用户主线程（通过 Run.log，已在 with_api 锁内串行化）
#   - next_system_step: Monitor Timer 线程（独立递增 system_step，与 next_step 互不干扰）
#   - define_scalar/define_media: Consumer 线程（单线程消费队列，串行）
# 如果未来引入多线程并发写入此对象，需要重新评估线程安全。
@dataclass
class RunMetrics:
    _global_step: int = 0
    _global_system_step: int = 0
    _metrics: Dict[str, Union[ScalarMetric, MediaMetric]] = field(default_factory=dict)

    def next_system_step(self) -> int:
        """
        获取下一个全局系统步数，用于系统内部监控指标
        """
        self._global_system_step += 1
        return self._global_system_step

    def next_step(self, user_step: Optional[int] = None) -> int:
        """
        获取下一个全局步数
        在设计上我们允许用户在log时乱序设置step，但是global_step永远是最大的或者自增的那个，
        因此我们需要一个方法来获取当前的global_step，并且保证global_step是自增的
        """
        if user_step is not None:
            self._global_step = max(self._global_step, user_step)
            return user_step
        self._global_step += 1
        return self._global_step

    def define_scalar(self, *, key: str, column: ColumnRecord) -> ScalarMetric:
        """
        定义一个标量指标
        :param key: 指标键
        :param column: 指标列记录
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        scalar_metric = ScalarMetric(_column=column)
        self._metrics[key] = scalar_metric
        return scalar_metric

    def define_media(self, *, key: str, column: ColumnRecord, path: Path) -> MediaMetric:
        """
        定义媒体指标
        :param key: 指标键
        :param column: 指标列记录
        :param path: 媒体存储路径，绝对路径
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        media_metric = MediaMetric(_column=column, path=path)
        self._metrics[key] = media_metric
        return media_metric

    def get(self, key: str) -> Optional[Union[ScalarMetric, MediaMetric]]:
        """
        根据key获取指标
        :param key: 指标键
        :return: 指标对象，如果不存在则返回None
        """
        return self._metrics.get(key)
