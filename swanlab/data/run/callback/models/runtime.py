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
from typing import Optional, Any
import json
import yaml
import os


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

    @abstractmethod
    def to_dict(self):
        """
        将信息转换为字典
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
        with open(os.path.join(path, self.name), "w", encoding="utf-8") as f:
            f.write(self.info)

    def dumps(self):
        return self.info

    def to_dict(self):
        raise NotImplementedError("RequirementInfo has no to_dict method")


class MetadataInfo(InfoWriter):
    """
    系统信息
    """

    def __init__(self, info: dict):
        super().__init__(info)
        self.name = "swanlab-metadata.json"

    def write(self, path: str):
        with open(os.path.join(path, self.name), "w", encoding="utf-8") as f:
            f.write(self.dumps())

    def dumps(self):
        return json.dumps(self.to_dict(), ensure_ascii=False)

    def to_dict(self):
        return self.info


class ConfigInfo(InfoWriter):
    """
    配置信息
    """

    def __init__(self, info: dict):
        super().__init__(info)
        self.name = "config.yaml"
        self.__data = None

    def write(self, path: str):
        with open(os.path.join(path, self.name), "w", encoding="utf-8") as f:
            f.write(self.dumps())

    def dumps(self):
        return yaml.dump(self.to_dict(), allow_unicode=True)

    def to_dict(self):
        """
        返回配置信息的字典形式
        """
        # 遍历每一个配置项，值改为value，增加desc和sort字段
        # 原因是序列化时可能会丢失key的排序信息，所以这里增加一个sort字段
        # 没有在__init__中直接修改是因为可能会有其他地方需要原始数据，并且会丢失一些性能
        if self.__data is not None:
            return self.__data
        self.__data = {
            k: {
                "value": v,
                "sort": i,
                "desc": ""
            } for i, (k, v) in enumerate(self.info.items())
        }
        return self.__data


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
        self.requirements: Optional[RequirementInfo] = RequirementInfo(
            requirements
        ) if requirements is not None else None

        self.metadata: Optional[MetadataInfo] = MetadataInfo(
            metadata
        ) if metadata is not None else None

        self.config: Optional[ConfigInfo] = ConfigInfo(
            config
        ) if config is not None else None
