#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 20:23
@File: callback.py
@IDE: pycharm
@Description:
    基本回调函数注册表，此时不考虑云端情况
"""
from typing import Union, Tuple
from swanlab.data.modules import DataType

NewKeyInfo = Union[None, Tuple[dict, Union[float, DataType], int, int]]
"""
新的key对象、数据类型、步数、行数
为None代表没有添加新的key
"""


class SwanLabRunCallback:

    def on_train_begin(self, *args, **kwargs):
        pass

    def on_metric_create(self, *args, **kwargs):
        pass

    def on_column_create(self, *args, **kwargs):
        pass

    def on_train_end(self, *args, **kwargs):
        pass
