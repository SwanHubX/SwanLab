"""
@author: cunyue
@file: metrics.py
@time: 2026/3/12 01:29
@description: SwanLab 运行时指标管理
"""

import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Literal, Optional, Union

from swanlab.sdk.internal.pkg import console

from .transformer import TransformType

# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


# 使用解包的方式传入参数
@dataclass(**DATACLASS_KWARGS)
class ScalarMetric:
    # 必须是普通实例变量，初始值为 None
    _type: TransformType
    # 指标对应的单实验图表索引
    _chart: str
    # 指标对应的图表名称，默认为列名称
    _chart_name: Optional[str] = None
    # 指标名称，默认为列名称
    _name: Optional[str] = None
    # 指标颜色
    _color: Optional[str] = None
    # x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
    _x_axis: Optional[Union[str, Literal["_step", "_relative_time"]]] = "_step"
    # 是否为系统指标
    _system: bool = False
    # 指标最新值
    latest: Optional[Union[float, int]] = None
    # 指标的最大值
    max: Optional[Union[float, int]] = None
    # 指标的最小值
    min: Optional[Union[float, int]] = None

    def update(self, _type: "TransformType", value: Union[float, int]) -> bool:
        """
        更新指标值
        :param _type: 传入的指标类型
        :param value: 指标值
        :return: 是否成功更新
        """
        # 首次写入逻辑：如果当前没有类型，则锁定为传入的类型
        if self._type is None:
            self._type = _type

        # 后续写入逻辑：校验类型是否一致
        elif self._type != _type:
            console.error(f"Invalid type for scalar metric: expected {self._type}, got {_type}")
            return False

        # 更新值
        self.latest = value
        if self.max is None or value > self.max:
            self.max = value
        if self.min is None or value < self.min:
            self.min = value
        return True


@dataclass(**DATACLASS_KWARGS)
class MediaMetric:
    # 媒体类型
    _type: TransformType
    # 媒体存储路径
    path: Path

    def update(self, _type: TransformType, path: Path) -> bool:
        """
        更新媒体信息
        :param _type: 媒体类型
        :param path: 媒体存储路径
        :return: 是否成功更新
        """
        if _type != self._type:
            console.error(f"Invalid type for media metric: expected {_type}, got {_type}")
            return False
        self.path = path
        return True


# 指标状态，实验运行过程中不断更新
@dataclass
class RunMetrics:
    _global_step: int = 0
    _step_lock: threading.Lock = field(default_factory=threading.Lock)
    _metrics: Dict[str, Union[ScalarMetric, MediaMetric]] = field(default_factory=dict)

    def next_step(self, user_step: Optional[int] = None) -> int:
        """
        线程安全地获取下一个 global_step。
        如果用户传入了 step，则同步更新内部的 global_step 最大值。
        """
        with self._step_lock:
            if user_step is not None:
                # 如果用户显式传了 step，如果step大于当前global_step，则更新global_step，否则不更新
                if user_step > self._global_step:
                    self._global_step = user_step
                return user_step
            else:
                # 隐式递增
                self._global_step += 1
                return self._global_step

    @property
    def global_step(self) -> int:
        """只读属性"""
        return self._global_step
