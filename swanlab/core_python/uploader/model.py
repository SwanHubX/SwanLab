"""
@author: cunyue
@file: model.py
@time: 2025/6/16 14:12
@description: 定义上传模型
"""

import json
from datetime import datetime
from enum import Enum
from typing import List, Optional, TypedDict, Literal

from swanlab.data.modules import MediaBuffer
from swanlab.toolkit import ColumnClass, ColumnConfig, LogContent


class ColumnModel:
    """
    列信息上传模型
    """

    def __init__(
        self,
        key,
        name: Optional[str],
        cls: ColumnClass,
        typ: str,
        config: Optional[ColumnConfig],
        section_name: Optional[str],
        section_type: Optional[str],
        error: dict = None,
    ):
        """
        key: 键
        name: 键的名称
        cls: 键的类别
        typ: 键的类型
        config: 键的配置
        section_name: 键所在的section的名称
        section_type: 键所在的section的类型
        error: 错误信息
        """
        self.key = key
        self.name = name
        self.cls = cls
        self.typ = typ
        self.config = config
        self.section_name = section_name
        self.section_type = section_type
        self.error = error

    def to_dict(self):
        """
        序列化为Dict,传递给后端
        """
        d = {
            "class": self.cls,
            "type": self.typ,
            "key": self.key,
            "name": self.name,
            "error": self.error,
            "sectionName": self.section_name,
            "sectionType": self.section_type,
        }
        if self.name is None:
            d.pop("name")
        if self.error is None:
            d.pop("error")
        if self.section_name is None:
            d.pop("sectionName")
        if self.section_type is None:
            d.pop("sectionType")
        if self.config is None:
            return d
        # 将额外的图表配置信息加入
        if self.config.y_range is not None:
            d["yRange"] = self.config.y_range
        if self.config.chart_name is not None:
            d["chartName"] = self.config.chart_name
        if self.config.chart_index is not None:
            d["chartIndex"] = self.config.chart_index
        if self.config.metric_name is not None:
            d["metricName"] = self.config.metric_name
        if self.config.metric_color is not None:
            d["metricColors"] = self.config.metric_color
        return d


class MetricType(Enum):
    """
    指标类型枚举
    """

    SCALAR = "scalar"
    """
    标量指标
    """
    MEDIA = "media"
    """
    媒体指标
    """
    LOG = "log"
    """
    日志数据
    """


class MediaModel:
    """
    媒体指标信息上传模型
    """

    type = MetricType.MEDIA

    def __init__(
        self,
        metric: dict,
        key: str,
        key_encoded: str,
        step: int,
        epoch: int,
        buffers: List[MediaBuffer] = None,
    ):

        # -------------------------- 🤡这里是一点小小的💩 --------------------------
        # 要求上传时的文件路径必须带key_encoded前缀
        if buffers is not None:
            metric = json.loads(json.dumps(metric))
            for i, d in enumerate(metric["data"]):
                metric["data"][i] = "{}/{}".format(key_encoded, d)
        # ------------------------------------------------------------------------

        self.metric = metric
        self.step = step
        self.epoch = epoch
        self.key = key
        """
        真实的指标名称
        """
        self.key_encoded = key_encoded
        """
        编码后路径安全的指标名称
        """
        self.buffers = buffers
        """
        原始数据，可能为None
        """

    def to_dict(self):
        """
        序列化
        """
        return {
            **self.metric,
            "key": self.key,
            "index": self.step,
            "epoch": self.epoch,
        }


class ScalarModel:
    """
    标量指标信息上传模型
    """

    type = MetricType.SCALAR

    def __init__(self, metric: dict, key: str, step: int, epoch: int):
        self.metric = metric
        self.key = key
        self.step = step
        self.epoch = epoch

    def to_dict(self):
        """
        序列化
        """
        return {
            **self.metric,
            "key": self.key,
            "index": self.step,
            "epoch": self.epoch,
        }


class FileModel:
    """
    运行时文件信息上传模型
    """

    def __init__(
        self,
        requirements: str = None,
        metadata: dict = None,
        config: dict = None,
        conda: str = None,
    ):
        self.requirements = requirements
        self.metadata = metadata
        self.config = config
        self.conda = conda
        # 主要用于去重，保留最新的文件
        self.create_time = datetime.now()

    @classmethod
    def create(cls, file_models: List["FileModel"]) -> "FileModel":
        """
        比较若干个FileModel，获取最新的FileModel，并且保证其内部属性不为None
        """
        # 按照时间排序，倒叙排列，这意味着最新的排在第一个
        file_models = sorted(file_models, key=lambda x: x.create_time, reverse=True)
        lr, lm, lc, lo = None, None, None, None
        for file_model in file_models:
            lr = file_model.requirements if lr is None else lr
            lm = file_model.metadata if lm is None else lm
            lc = file_model.config if lc is None else lc
            lo = file_model.conda if lo is None else lo
            if lr is not None and lm is not None and lc is not None and lo is not None:
                break
        return FileModel(lr, lm, lc, lo)

    def to_dict(self):
        """
        序列化，会删除为None的字段
        """
        d = {
            "requirements": self.requirements,
            "metadata": self.metadata,
            "config": self.config,
            "conda": self.conda,
        }
        return {k: v for k, v in d.items() if v is not None}

    @property
    def empty(self):
        """
        是否为空
        """
        return all(getattr(self, attr) is None for attr in ('requirements', 'metadata', 'config', 'conda'))


class LogModel(TypedDict):
    """
    日志信息上传模型
    """

    level: Literal['INFO', 'WARN', 'ERROR']
    """
    日志级别
    """
    contents: List[LogContent]
    """
    当前日志级别的日志内容
    """


__all__ = [
    "LogModel",
    "ColumnModel",
    "ScalarModel",
    "MediaModel",
    "FileModel",
]
