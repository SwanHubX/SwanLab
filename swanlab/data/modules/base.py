#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-28 22:12:37
@File: swanlab/database/modules/base.py
@IDE: vscode
@Description:
    SwanLab数据基类
"""
from abc import ABC, abstractmethod
from ..settings import SwanDataSettings
from .chart import Chart
from urllib.parse import quote


class BaseType(ABC):
    """SwanLab数据基类，完成数据解析，拿到图表类型、namespace和解析后的数据"""

    def __init__(self, value, *args, **kwargs):
        # 保存的key
        self.__tag: str = None
        # 保存的step
        self.__step: int = None
        # 受支持的图表类型，不可修改
        self.__chart: Chart = Chart()
        # 当前运行时配置
        self.__settings: SwanDataSettings = None
        # 保存的原始值
        self.value = value
        # 标志是否已经被解构
        self.__extracted: bool = False
        # get_data方法是否被执行
        self.__is_get_data: bool = False
        # get_data的返回值
        self.__data = None

    def __iter__(self):
        # 返回一个迭代器，通常是自身
        if self.__extracted:
            raise StopIteration("This object has been extracted")
        return self

    def __next__(self):
        """解构对象，用法是 namespace, chart_type, config = object"""
        return self.get_namespace(), self.get_chart_type(), "step", self.get_config()

    @property
    def chart(self) -> Chart:
        return self.__chart

    @property
    def convert(self):
        if self.__is_get_data is False:
            self.__data = self.get_data()
            self.__is_get_data = True
        return self.__data

    @abstractmethod
    def get_data(self, *args, **kwargs):
        """获取data数据，将返回一个可被JSON序列化的数据"""
        pass

    def get_namespace(self, *args, **kwargs) -> str:
        """获取namespace数据，将返回一个可被JSON序列化的数据"""
        return "default"

    @abstractmethod
    def get_chart_type(self, *args, **kwargs) -> str:
        """获取chart_type数据，应该返回一个字符串"""
        pass

    def get_config(self, *args, **kwargs) -> dict:
        """获取图表的config数据，应该返回一个字典"""
        return {}

    def get_more(self, *args, **kwargs):
        """代表当前步骤的此数据支持标注的更多内容，应该返回一个字典，或者为None"""
        return None

    @abstractmethod
    def expect_types(self, *args, **kwargs) -> list:
        """输入期望的数据类型，list类型，例如["int", "float"]"""
        pass

    def get_data_list(self):
        """单step多数据时，返回多个数据的列表"""
        datas = []
        for i in self.value:
            try:
                i.settings = self.settings
            except AttributeError:
                pass
            try:
                i.tag = self.tag
            except AttributeError:
                pass
            try:
                i.step = self.step
            except AttributeError:
                pass
            datas.append(i.get_data())
        return [i.get_data() for i in self.value]

    def get_more_list(self):
        """单step多数据时，返回多个数据更多信息"""
        return [i.get_more() for i in self.value]

    # ---------------------------------- 向子类暴露tag、step、create_time,他们的值只能修改一次 ----------------------------------
    @property
    def tag(self):
        return self.__tag

    @property
    def step(self):
        return self.__step

    @property
    def settings(self):
        return self.__settings

    @settings.setter
    def settings(self, value):
        if self.__settings is None:
            self.__settings = value
        else:
            raise AttributeError("settings can only be set once")

    @tag.setter
    def tag(self, value):
        if self.__tag is None:
            # 将tag进行url编码,但是可能已经被编码了，所以需要设置safe="%2F"
            self.__tag = quote(value, safe="%2F")
        else:
            raise AttributeError("tag can only be set once")

    @step.setter
    def step(self, value):
        if self.__step is None:
            self.__step = value
        else:
            raise AttributeError("step can only be set once")
