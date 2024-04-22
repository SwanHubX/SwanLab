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

    def __init__(self, key, column_type: str, error: dict = None):
        """
        :param key: 列名称
        :param column_type: 列类型，'FLOAT', 'IMAGE', 'AUDIO', 'TEXT', 'ANY'，必须为大写，如果传入 'DEFAULT'，则会转为 'FLOAT'
        :param error: 错误信息，如果错误信息不为None，column_type必须为'ANY'
        """
        self.key = key
        if column_type == "DEFAULT":
            column_type = "FLOAT"
        self.column_type = column_type
        self.error = error

    def to_dict(self):
        if self.error is not None:
            return {
                "key": self.key,
                "type": 'ANY',
                "error": self.error
            }
        else:
            return {
                "key": self.key,
                "type": self.column_type,
            }
