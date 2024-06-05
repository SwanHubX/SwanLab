#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/5 16:55
@File: runtime.py
@IDE: pycharm
@Description:
    运行时信息模型
"""
from abc import ABC, abstractmethod
import json
import yaml


class InfoWriter(ABC):

    def __init__(self, info):
        self.info = info

    @abstractmethod
    def write(self, path: str):
        """
        写入到本地文件的方法
        :param path: 写入的路径，不包含文件名
        """
        pass

    @abstractmethod
    def dumps(self):
        """
        将信息转换为字符串
        """
        pass


class RequirementInfo(InfoWriter):
    """
    依赖信息
    """

    def __init__(self, info: str):
        super().__init__(info)
        self.name = "requirements.txt"

    def write(self, path: str):
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.info)

    def dumps(self):
        return self.info


class MetadataInfo(InfoWriter):
    """
    系统信息
    """

    def __init__(self, info: dict):
        super().__init__(info)
        self.name = "swanlab-metadata.json"

    def write(self, path: str):
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.dumps())

    def dumps(self):
        return json.dumps(self.info, ensure_ascii=False)


class ConfigInfo(InfoWriter):
    """
    配置信息
    """

    def __init__(self, info: dict):
        super().__init__(info)
        self.name = "config.yaml"

    def write(self, path: str):
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.dumps())

    def dumps(self):
        return yaml.dump(self.info, allow_unicode=True)


class RuntimeInfo:
    """
    运行时信息，包括系统信息，依赖信息等
    如果某些信息为None，代表没有传入配置
    """

    def __init__(self, requirements: str = None, metadata: dict = None, config: dict = None):
        """
        :param requirements: python依赖信息
        :param metadata: 系统信息
        :param config: 上传的配置信息
        """
        self.requirements = RequirementInfo(requirements) if requirements is not None else None
        self.metadata = MetadataInfo(metadata) if metadata is not None else None
        self.config = ConfigInfo(config) if config is not None else None
