#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/22 16:42
@File: model.py
@IDE: pycharm
@Description:
    上传请求模型
"""
from enum import Enum
from typing import List
from swanlab.data.modules import MediaBuffer
from datetime import datetime


class ColumnModel:
    """
    列信息上传模型
    """

    def __init__(self, key, column_type: str, error: dict = None):
        """
        :param key: 列名称
        :param column_type: 列类型，'FLOAT', 'IMAGE', 'AUDIO', 'TEXT'，必须为大写，如果传入 'DEFAULT'，则会转为 'FLOAT'
        :param error: 错误信息，如果错误信息不为None
        """
        self.key = key
        self.column_type = column_type
        self.error = error

    def to_dict(self):
        """
        序列化为Dict
        """
        return {
            "key": self.key,
            "type": self.column_type,
        } if self.error is None else {
            "key": self.key,
            "type": self.column_type,
            "error": self.error
        }


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
        buffers: List[MediaBuffer] = None
    ):
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
            "epoch": self.epoch
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
            "epoch": self.epoch
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
        create_time: datetime = None
    ):
        self.requirements = requirements
        self.metadata = metadata
        self.config = config
        self.create_time = datetime.now()
        """
        主要用于去重，保留最新的文件
        """

    @classmethod
    def create(cls, file_models: List["FileModel"]) -> "FileModel":
        """
        比较若干个FileModel，获取最新的FileModel，并且保证其内部属性不为None
        """
        # 按照时间排序，倒叙排列，这意味着最新的排在第一个
        file_models = sorted(file_models, key=lambda x: x.create_time, reverse=True)
        lr, lm, lc = None, None, None
        for file_model in file_models:
            lr = file_model.requirements if lr is None else lr
            lm = file_model.metadata if lm is None else lm
            lc = file_model.config if lc is None else lc
            if lr is not None and lm is not None and lc is not None:
                break
        return FileModel(lr, lm, lc)

    def to_dict(self):
        """
        序列化，会删除为None的字段
        """
        d = {"requirements": self.requirements, "metadata": self.metadata, "config": self.config}
        return {k: v for k, v in d.items() if v is not None}
