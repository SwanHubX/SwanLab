#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/22 16:42
@File: model.py
@IDE: pycharm
@Description:
    上传请求模型
"""


class ColumnModel:
    """
    列信息上传模型
    """

    def __init__(self, key, column_type, error: dict = None):
        """
        :param key: 列名称
        :param column_type: 列类型，'FLOAT', 'IMAGE', 'AUDIO', 'TEXT', 'ANY'
        :param error: 错误信息，如果错误信息不为None，column_type必须为'ANY'
        """
        self.key = key
        self.column_type = column_type
        self.error = error

    def to_dict(self):
        return {
            "key": self.key,
            "columnType": self.column_type,
            "error": self.error
        }
