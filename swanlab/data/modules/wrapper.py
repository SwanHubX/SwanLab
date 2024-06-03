#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:39
@File: wrapper.py
@IDE: pycharm
@Description:
    包装器
"""
from .base import BaseType, MediaType, ParseResult
from typing import Union, List, Optional
from swanlab.error import DataTypeError
from .line import Line


class ErrorInfo:
    """
    DataWrapper转换时的错误信息
    """

    def __init__(self, expected: str, got: str, chart: BaseType.Chart):
        """
        :param expected: 期望的数据类型
        :param got: 实际的数据类型
        :param chart: 当前错误数据对应的图表类型
        """
        self.expected = expected
        self.got = got
        self.chart = chart

    def dict(self):
        return {"expected": self.expected, "got": self.got}


class DataWrapper:
    """
    数据包装器，负责解析数据，返回相应的数据类型

    一个数据包装器对应swanlab.log字典内的一个key的值
    """

    def __init__(self, key: str, data: Union[Line, List[MediaType]]):
        """
        :param key: key名称，为编码前的名称
        :param data: log的数据
        """
        self.__key = key
        self.__data = [data] if isinstance(data, Line) else data
        self.__error = None
        self.__result: Optional[ParseResult] = None
        # 保证data内所有数据类型相同，否则返回TypeError
        t = self.__data[0].__class__
        for i in self.__data:
            if i.__class__ != t:
                raise TypeError(f"All data in DataWrapper must be the same type, got {t} and {i.__class__}")

    @property
    def is_line(self) -> bool:
        """
        是否是Line类型
        """
        return isinstance(self.__data[0], Line) and len(self.__data) == 1

    @property
    def type(self):
        """
        数据类型
        """
        return self.__data[0].__class__

    @property
    def parsed(self) -> bool:
        """
        是否已经解析过数据
        """
        return self.__result is not None or self.__error is not None

    @property
    def error(self) -> Optional[ErrorInfo]:
        """
        解析时候的错误信息
        """
        return self.__error

    def parse(self, *args, **kwargs) -> Optional[ParseResult]:
        """
        将数据解析成对应的数据类型
        *args, **kwargs: 在解析数据类型之前，需要注入的信息
        """
        if self.parsed:
            return self.__result
        result = ParseResult()
        # ---------------------------------- 特别处理Line类型 ----------------------------------

        # [Line]
        if isinstance(self.type, Line):
            d = self.__data[0]
            result.section = d.get_section()
            result.chart = d.get_chart()
            result.write_handler = d.get_raw_write_handler()
            if len(self.__data) > 1:
                self.__error = ErrorInfo("Line", "list(Line)", result.chart)
            else:
                d.inject(*args, **kwargs)
                try:
                    result.data, _ = d.parse()
                except DataTypeError as e:
                    self.__error = ErrorInfo(e.expected, e.got, result.chart)
            self.__result = result
            return self.__result

        # ---------------------------------- 处理其他类型 ----------------------------------
        # [MediaType]
        data, raw, result.more, result.config = [], [], [], []
        for i in self.__data:
            try:
                i.inject(*args, **kwargs)
                d, r = i.parse()
            except DataTypeError as e:
                self.__error = ErrorInfo(e.expected, e.got, result.chart)
                return None
            data.append(d)
            raw.append(r)
            result.more.append(i.get_more())
            result.config.append(i.get_config())
        self.__result = result
        return self.__result
