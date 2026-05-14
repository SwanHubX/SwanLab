"""
@author: cunyue
@file: __init__.py
@time: 2026/5/5 18:53
@description: SwanLab指标管理模块，负责在本地管理指标上下文
指标上下文统一由Core处理的好处是不需要SDK（多语言）前端重复管理指标上下文，前端只管生产，后端负责拦截和管理
"""

import math
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Set, Tuple, Union, cast

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg import builder
from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.typings.core_python.api.experiment import ResumeExperimentSummaryType

__all__ = ["RunMetrics"]

# 开启 slots 可以减少内存占用，提高性能，但是仅适用于 Python 3.10 及以上版本
# 因此我们动态判断版本：如果是 3.10 及以上，开启 slots；否则传空字典
DATACLASS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


@dataclass(**DATACLASS_KWARGS)
class BaseMetric(ABC):
    _column: ColumnRecord
    min_step: int
    steps: Set[int]

    @property
    def column_type(self) -> ColumnType:
        return self._column.column_type

    def ensure_type_match(self, metric_type: ColumnType):
        if self.column_type != metric_type:
            raise TypeError(
                f"Metric '{self._column.column_key}' has already been defined as "
                f"{ColumnType.Name(self.column_type)}, not {ColumnType.Name(metric_type)}."
            )

    def try_accept_step(self, step: int) -> bool:
        """
        检查当前 step 是否应跳过；如果是新的合法 step，则登记下来。

        当出现以下情况时返回 False：
        - step 小于等于 min_step
        - step 已经记录过

        当 step 可接受且此前未记录时：
        - 将 step 加入已记录集合
        - 返回 True

        :param step: 当前步数
        :return: 是否应跳过该 step
        """
        if step <= self.min_step:
            return False

        # 取消下面的处理是实现 https://github.com/SwanHubX/SwanLab/issues/1576 的前置条件，我们需要等待后端准备好
        if step in self.steps:
            return False
        self.steps.add(step)
        return True

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

    def update(self, data_record: Any):
        """
        更新标量指标值，为了兼容父类的类型定义，这里使用 Any 类型，并用cast 转换为 ScalarRecord 类型
        :param data_record: 标量指标记录
        """
        value = cast(ScalarRecord, data_record).value.number
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

    def ensure_type_match(self, metric_type: ColumnType):
        """
        Resume 模式下我们可能暂时无法知道媒体指标的具体类型，因此在 ensure_type_match 中直接 pass
        """
        if self.column_type == ColumnType.COLUMN_TYPE_UNSPECIFIED:
            console.debug(
                f"Media metric '{self._column.column_key}' has not been defined, maybe it's defined in resume"
            )
            return
        BaseMetric.ensure_type_match(self, metric_type)

    def update(self, data_record: Any):
        """
        更新媒体指标值，媒体指标目前没什么可更新的，因此直接pass
        """
        pass


# 指标状态，实验运行过程中不断更新
# 注意：此类不加锁保护，如果未来引入多线程并发写入此对象，需要重新评估线程安全。
@dataclass
class RunMetrics:
    _metrics: Dict[str, Union[ScalarMetric, MediaMetric]] = field(default_factory=dict)

    @classmethod
    def new(cls, data: Optional[ResumeExperimentSummaryType], ctx: CoreContext) -> Tuple["RunMetrics", int, int, int]:
        """
        根据后端返回的实验摘要，生成 RunMetrics 对象
        :param data: 后端返回的实验摘要
        :param ctx: 运行上下文
        :return: RunMetrics 对象, 起始步数
        """
        console_epoch, global_step, global_system_step = 0, -1, -1
        if data is None:
            return cls(), console_epoch, global_step, global_system_step
        # 全局起始步数，系统起始步数，分别用于用户侧起始步数和系统采集侧起始步数
        # 1. 获取起始的终端记录行数
        if "log" in data and data["log"] is not None and len(data["log"]) > 0:
            console_epoch = data["log"][0]["step"]
        metrics = cls()
        # 2. 获取媒体指标记录
        if "media" in data and data["media"] is not None:
            for media in data["media"]:
                media_key, media_step = media["key"], media["step"]
                if media_step > global_step:
                    global_step = media_step
                column_record = builder.build_resume_column(media_key, media=True)
                path = ctx.media_dir / "unknown"
                metrics.define_media(key=media_key, column=column_record, path=path, min_step=media_step)
        # 3. 获取标量指标记录
        if "scalar" in data and data["scalar"] is not None:
            for scalar in data["scalar"]:
                scalar_key, scalar_step = scalar["key"], scalar["step"]
                # 系统列和自定义列有不同的起始步数
                is_system = helper.is_system_key(scalar_key)
                if is_system:
                    global_system_step = max(scalar_step, global_system_step)
                else:
                    global_step = max(scalar_step, global_step)
                column_record = builder.build_resume_column(scalar_key, system=is_system)
                metrics.define_scalar(key=scalar_key, column=column_record, min_step=scalar_step)
        return metrics, console_epoch, global_step, global_system_step

    def define_scalar(self, *, key: str, column: ColumnRecord, min_step: int = -1) -> ScalarMetric:
        """
        定义一个标量指标
        :param key: 指标键
        :param column: 指标列记录
        :param min_step: 最小步数，代表用户无法再写入此步数之前的数据，这在阻止用户写入step小于0、resume时拒绝一定大小的step十分有用
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        scalar_metric = ScalarMetric(_column=column, steps=set(), min_step=min_step)
        self._metrics[key] = scalar_metric
        return scalar_metric

    def define_media(self, *, key: str, column: ColumnRecord, path: Path, min_step: int = -1) -> MediaMetric:
        """
        定义媒体指标
        :param key: 指标键
        :param column: 指标列记录
        :param path: 媒体存储路径，绝对路径
        :param min_step: 最小步数，代表用户无法再写入此步数之前的数据，这在阻止用户写入step小于0、resume时拒绝一定大小的step十分有用
        """
        assert key not in self._metrics, f"Metric '{key}' already exists."
        media_metric = MediaMetric(_column=column, path=path, steps=set(), min_step=min_step)
        self._metrics[key] = media_metric
        return media_metric

    def get(self, key: str) -> Optional[Union[ScalarMetric, MediaMetric]]:
        """
        根据key获取指标
        :param key: 指标键
        :return: 指标对象，如果不存在则返回None
        """
        return self._metrics.get(key)
