#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 02:32
@File: __init__.py.py
@IDE: pycharm
@Description:
    浮点数模块
    解析浮点数
    1. 正常情况下返回浮点数, None
    2. 如果传入的数为不可转变为浮点数的类型，抛出DataTypeError异常
    3. 如果传入NaN，返回NaN字符串，None
    4. 如果传入Infinity，返回INF字符串，None
"""
from swanlab.error import DataTypeError
from typing import Protocol, runtime_checkable
from swankit.core import BaseType, DataSuite as D


@runtime_checkable
class FloatConvertible(Protocol):
    """
    实现了__float__方法的类
    """

    def __float__(self) -> float: ...


class Line(BaseType):
    nan = "NaN"
    inf = "INF"

    def __init__(self, value):
        super().__init__()
        self.value = value

    def parse(self):
        # 如果是nan
        try:
            t = float(self.value)
            if D.is_nan(t):
                return self.nan, None
            # 如果是inf
            if D.is_inf(t):
                return self.inf, None
            return t, None
        except (ValueError, TypeError):
            raise DataTypeError('float', type(self.value).__name__)

    def get_chart(self):
        return self.Chart.LINE
