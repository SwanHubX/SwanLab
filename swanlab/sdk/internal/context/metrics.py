"""
@author: cunyue
@file: metrics.py
@time: 2026/3/12 01:29
@description: SwanLab 运行时指标管理
"""

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Union

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.sdk.internal.pkg import console

# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


# 使用解包的方式传入参数
@dataclass(**DATACLASS_KWARGS)
class ScalarMetric:
    # 指标列记录
    _column: ColumnRecord
    # 指标最新值
    latest: Optional[Union[float, int]] = None
    # 指标的最大值
    max: Optional[Union[float, int]] = None
    # 指标的最小值
    min: Optional[Union[float, int]] = None

    def update(self, value: Union[float, int]) -> bool:
        """
        更新标量指标值
        :param value: 标量指标值
        :return: 是否成功更新
        """
        self.latest = value
        if self.max is None or value > self.max:
            self.max = value
        if self.min is None or value < self.min:
            self.min = value
        return True

    @property
    def type(self):
        return self._column.column_type


@dataclass(**DATACLASS_KWARGS)
class MediaMetric:
    # 指标列记录
    _column: ColumnRecord
    # 媒体存储路径，绝对路径
    path: Path

    @property
    def type(self):
        return self._column.column_type


# 指标状态，实验运行过程中不断更新
# 注意：此类不加锁保护。设计上对 RunMetrics 的写操作发生在以下线程：
#   - next_step: 用户主线程（通过 Run.log，已在 with_api 锁内串行化）
#   - next_system_step: Monitor Timer 线程（独立递增 system_step，与 next_step 互不干扰）
#   - update_scalar: Consumer 线程（单线程消费队列，与 define_scalar/define_media 串行）
#   - define_scalar/define_media: Consumer 线程（同上）
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

    def update_scalar(self, key: str, value: Union[float, int]):
        """
        更新标量指标状态
        :param key: 指标键
        :param value: 标量值
        """
        scalar = self._metrics.get(key)
        assert scalar is not None, f"Metric '{key}' does not exist."
        assert isinstance(scalar, ScalarMetric), f"Metric '{key}' is not a scalar metric."
        if math.isnan(value) or math.isinf(value):
            console.debug(f"Invalid scalar value: {value} for metric '{key}', ignored when updating.")
            return
        scalar.latest = value
        if scalar.max is None or value > scalar.max:
            scalar.max = value
        if scalar.min is None or value < scalar.min:
            scalar.min = value

    def define_scalar(self, key: str, column: ColumnRecord):
        """
        定义一个标量指标
        :param key: 指标键
        :param column: 指标列记录
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        self._metrics[key] = ScalarMetric(_column=column)

    def define_media(self, key: str, column: ColumnRecord, path: Path):
        """
        定义媒体指标
        :param key: 指标键
        :param column: 指标列记录
        :param path: 媒体存储路径，绝对路径
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        self._metrics[key] = MediaMetric(_column=column, path=path)

    def ensure_defined_as(self, key: str, metric_type: ColumnType) -> bool:
        """
        判断指标是否已定义，并确保其类型与预期一致。

        :param key: 指标键
        :param metric_type: 期望的指标类型
        :return:
            - False: 指标尚未定义
            - True: 指标已定义且类型一致
        :raise TypeError: 指标已定义但类型不匹配
        """
        if key not in self._metrics:
            return False

        metric = self._metrics[key]
        if metric.type != metric_type:
            raise TypeError(
                f"Metric '{key}' has already been defined as "
                f"{ColumnType.Name(metric.type)}, not {ColumnType.Name(metric_type)}."
            )
        return True
