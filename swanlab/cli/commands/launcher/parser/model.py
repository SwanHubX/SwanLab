#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:52
@File: model.py
@IDE: pycharm
@Description:
    基础模型
"""
import click
from typing import Type, Any
from abc import ABC, abstractmethod
import os


class LaunchParser(ABC):

    def __init__(self, config: dict, path: str):
        self.config = config
        """
        配置文件内容
        """
        self.dirpath = os.path.dirname(path)
        """
        配置文件所在目录
        """

    @staticmethod
    def should_be(n: str, v: Any, t: Type, none: bool = False) -> Any:
        """
        检查参数是否符合预期
        :param n: 参数名
        :param v: 参数值
        :param t: 预期类型
        :param none: 是否允许为None
        """
        if none and v is None:
            return None
        if v is None and not none:
            raise click.BadParameter(f'{n} should not be None')
        if not isinstance(v, t):
            raise click.BadParameter(f'{n} should be {t}, not {type(v)}')
        return v

    def should_file_exist(self, n: str, p: str):
        """
        检查文件是否存在，必须是文件
        :param n: 参数名
        :param p: 文件路径
        """
        p = os.path.join(self.dirpath, p)
        if not os.path.exists(p):
            raise click.FileError(p, hint=f'{n} not found: {p}')
        if not os.path.isfile(p):
            raise click.FileError(p, hint=f'{n} should be a file: {p}')

    @staticmethod
    def should_in_values(n: str, v: Any, ls: list) -> Any:
        """
        检查参数是否在指定范围内
        :param n: 参数名
        :param v: 参数值
        :param ls: 范围
        """
        if v not in ls:
            raise click.BadParameter(f'{n} should be in {ls}, not {v}')
        return v

    @staticmethod
    def should_equal_keys(n: str, v: dict, keys: list) -> Any:
        """
        确保字典的key在指定范围内
        :param n: 参数名
        :param v: 参数值
        :param keys: 范围
        """
        for k in v.keys():
            if k not in keys:
                raise click.BadParameter(f'Unknown key: {n}.{k}')
        return v

    @classmethod
    @abstractmethod
    def __type__(cls):
        """
        返回当前解析器的名称，对应到kind参数
        """
        pass

    @abstractmethod
    def __dict__(self) -> dict:
        """
        最终向后端发布任务的data数据
        """
        pass

    @abstractmethod
    def parse_spec(self, spec: dict):
        """
        解析spec数据
        """
        pass

    @abstractmethod
    def parse_metadata(self, metadata: dict):
        """
        解析metadata数据
        """
        pass

    def parse(self):
        """
        解析整个配置文件
        """
        metadata = self.should_be('metadata', self.config.get('metadata'), dict)
        self.should_equal_keys('metadata', metadata, ['name', 'desc', 'combo'])
        self.parse_metadata(metadata)
        spec = self.should_be('spec', self.config.get('spec'), dict)
        self.should_equal_keys('spec', spec, ['python', 'entry', 'volumes', 'exclude'])
        self.parse_spec(spec)

    @abstractmethod
    def run(self):
        """
        执行具体的操作
        """
        pass

    @abstractmethod
    def dry_run(self):
        """
        单纯向用户展示即将执行的操作，而不实际应用
        """
        pass
