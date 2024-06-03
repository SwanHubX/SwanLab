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

    def __init__(self, metric: dict, key: str, key_encoded: str, step: int, epoch: int, raw: bytes = None):
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
        self.raw = raw
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
