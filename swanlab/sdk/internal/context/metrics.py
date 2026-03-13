"""
@author: cunyue
@file: metrics.py
@time: 2026/3/12 01:29
@description: SwanLab 运行时指标管理
"""

import math
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Literal, Optional, Union

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.data import MediaTransferType, ScalarXAxisType

# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


# 使用解包的方式传入参数
@dataclass(**DATACLASS_KWARGS)
class ScalarMetric:
    # 指标对应的单实验图表索引
    _chart: Optional[str] = None
    # 指标对应的图表名称，默认为列名称
    _chart_name: Optional[str] = None
    # 指标名称，默认为列名称
    _name: Optional[str] = None
    # 指标颜色
    _color: Optional[str] = None
    # x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
    _x_axis: Union[str, Literal["_step", "_relative_time"]] = "_step"
    # 是否为系统指标
    _system: bool = False
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


@dataclass(**DATACLASS_KWARGS)
class MediaMetric:
    # 媒体类型
    _type: MediaTransferType
    # 媒体存储路径
    path: Path


# 指标状态，实验运行过程中不断更新
@dataclass
class RunMetrics:
    _global_step: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock)
    _metrics: Dict[str, Union[ScalarMetric, MediaMetric]] = field(default_factory=dict)

    def next_step(self, user_step: Optional[int] = None) -> int:
        """
        获取下一个全局步数
        """
        with self._lock:
            if user_step is not None:
                self._global_step = max(self._global_step, user_step)
                return user_step
            self._global_step += 1
            return self._global_step

    def has_metric(self, key: str) -> bool:
        """
        检查是否存在指定的指标
        :param key: 指标键
        :return: 是否存在该指标
        """
        with self._lock:
            return key in self._metrics

    def update_scalar(self, key: str, value: Union[float, int]):
        """
        更新标量指标状态
        :param key: 指标键
        :param value: 标量值
        """
        with self._lock:
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

    def define_scalar(
        self,
        key: str,
        name: Optional[str] = None,
        chart: Optional[str] = None,
        chart_name: Optional[str] = None,
        system: bool = False,
        color: Optional[str] = None,
        x_axis: Optional[ScalarXAxisType] = None,
    ):
        """
        定义一个标量指标
        :param key: 指标键
        :param name: 指标名称，默认为列名称
        :param chart: 单实验图表索引
        :param chart_name: 图表名称，默认为列名称
        :param system: 是否为系统指标
        :param color: 指标颜色
        :param x_axis: x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
        :return:
        """
        with self._lock:
            assert key not in self._metrics, f"Metric '{key}' already exists."
            x_axis = x_axis or "_step"
            self._metrics[key] = ScalarMetric(
                _chart=chart,
                _chart_name=chart_name,
                _name=name,
                _system=system,
                _color=color,
                _x_axis=x_axis,
            )

    def define_media(self, key: str, media_type: MediaTransferType, path: Path):
        """
        定义媒体指标
        :param key: 指标键
        :param media_type: 媒体类型
        :param path: 媒体存储路径
        """
        with self._lock:
            assert key not in self._metrics, f"Metric '{key}' already exists."
            self._metrics[key] = MediaMetric(_type=media_type, path=path)
