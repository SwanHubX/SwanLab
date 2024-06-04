#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:40
@File: config.py
@IDE: pycharm
@Description:
    SwanLabConfig 配置类
"""
from typing import Any
from collections.abc import Mapping
import yaml
import argparse
from ..settings import SwanDataSettings
from swanlab.log import swanlog
import datetime
import math


def json_serializable(obj: dict):
    """
    将传入的字典转换为JSON可序列化格式。
    """
    # 如果对象是基本类型，则直接返回
    if isinstance(obj, (int, float, str, bool, type(None))):
        if isinstance(obj, float) and math.isnan(obj):
            return "nan"
        if isinstance(obj, float) and math.isinf(obj):
            return "inf"
        return obj

    # 将日期和时间转换为字符串
    elif isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()

    # 对于列表和元组，递归调用此函数
    elif isinstance(obj, (list, tuple)):
        return [json_serializable(item) for item in obj]

    # 对于字典，递归调用此函数处理值，并将key转换为字典
    elif isinstance(obj, dict):
        return {str(key): json_serializable(value) for key, value in obj.items()}

    else:
        # 对于其他不可序列化的类型，转换为字符串表示
        return str(obj)


def thirdparty_config_process(data) -> dict:
    """
    对于一些特殊的第三方库的处理，例如omegaconf
    """
    # 如果是omegaconf的DictConfig，则转换为字典
    try:
        import omegaconf

        if isinstance(data, omegaconf.DictConfig):
            return omegaconf.OmegaConf.to_container(data, resolve=True, throw_on_missing=True)
    except Exception as e:
        pass

    # 如果是argparse的Namespace，则转换为字典
    if isinstance(data, argparse.Namespace):
        return vars(data)

    return data


class SwanLabConfig(Mapping):
    """
    The SwanConfig class is used for realize the invocation method of `run.config.lr`.
    """

    # 配置字典
    __config = dict()

    # 运行时设置
    __settings = dict()

    @property
    def _inited(self):
        return self.__settings.get("save_path") is not None

    def __init__(self, config: dict = None, settings: SwanDataSettings = None):
        """
        实例化配置类，如果settings不为None，说明是通过swanlab.init调用的，否则是通过swanlab.config调用的

        Parameters
        ----------
        settings : SwanDataSettings, optional
            运行时设置
        """
        if config is None:
            config = {}

        self.__config.update(config)
        self.__settings["save_path"] = settings.config_path if settings is not None else None
        self.__settings["should_save"] = settings.should_save if settings is not None else False
        if self._inited:
            self.save()

    @property
    def should_shave(self):
        return self.__settings.get("should_save")

    @staticmethod
    def __check_config(config: dict) -> dict:
        """
        检查配置是否合法，确保它可以被 JSON/YAML 序列化，并返回转换后的配置字典。
        """
        if config is None:
            return {}
        # config必须可以被json序列化
        try:
            # 第三方配置类型判断与转换
            config = thirdparty_config_process(config)

            # 将config转换为可被json序列化的字典
            config = json_serializable(dict(config))

            # 尝试序列化，如果还是失败就退出
            yaml.dump(config)

        except:
            raise TypeError(f"config: {config} is not a valid dict, which can be json serialized")

        return config

    @staticmethod
    def __check_private(name: str):
        """
        检查属性名是否是私有属性,如果是私有属性，抛出异常

        Parameters
        ----------
        name : str
            属性名

        Raises
        ----------
        AttributeError
            如果属性名是私有属性，抛出异常
        """
        methods = ["set", "get", "pop"]
        swanlog.debug(f"Check private attribute: {name}")
        if name.startswith("__") or name.startswith("_SwanLabConfig__") or name in methods:
            raise AttributeError("You can not get private attribute")

    def __setattr__(self, name: str, value: Any) -> None:
        """
        自定义属性设置方法。如果属性名不是私有属性，则同时更新配置字典并保存。
        允许通过点号方式设置属性，但不允许设置私有属性：
        ```python
        run.config.lr = 0.01  # 允许
        run.config._lr = 0.01 # 允许
        run.config.__lr = 0.01 # 不允许
        ```

        值得注意的是类属性的设置不会触发此方法
        """

        # 判断是否是私有属性
        name = str(name)
        self.__check_private(name)
        # 设置属性，并判断是否已经初始化，如果是，则调用保存方法
        self.__dict__[name] = value
        # 同步到配置字典
        self.__config[name] = value
        self.save()

    def __setitem__(self, name: str, value: Any) -> None:
        """
        以字典方式设置配置项的值，并保存，但不允许设置私有属性：
        ```python
        run.config["lr"] = 0.01  # 允许
        run.config["_lr"] = 0.01 # 允许
        run.config["__lr"] = 0.01 # 不允许
        ```
        """
        # 判断是否是私有属性
        name = str(name)
        self.__check_private(name)
        self.__config[name] = value
        self.save()

    def set(self, name: str, value: Any) -> None:
        """
        Explicitly set the value of a configuration item and save it. For example:

        ```python
        run.config.set("lr", 0.01)   # Allowed
        run.config.set("_lr", 0.01)  # Allowed
        run.config.set("__lr", 0.01) # Not allowed
        ```

        Parameters
        ----------
        name: str
            Name of the configuration item
        value: Any
            Value of the configuration item

        Raises
        ----------
        AttributeError
            If the attribute name is private, an exception is raised
        """
        name = str(name)
        self.__check_private(name)
        self.__config[name] = value
        self.save()

    def pop(self, name: str) -> bool:
        """
        Delete a configuration item; if the item does not exist, skip.

        Parameters
        ----------
        name : str
            Name of the configuration item

        Returns
        ----------
        bool
            True if deletion is successful, False otherwise
        """
        try:
            del self.__config[name]
            self.save()
            return True
        except KeyError:
            return False

    def get(self, name: str):
        """
        Get the value of a configuration item. If the item does not exist, raise AttributeError.

        Parameters
        ----------
        name : str
            Name of the configuration item

        Returns
        ----------
        value : Any
            Value of the configuration item

        Raises
        ----------
        AttributeError
            If the configuration item does not exist, an AttributeError is raised
        """
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not retrieved '{name}' in the config of the current experiment")

    def clean(self):
        """
        清空配置字典
        """
        self.__config.clear()

    def update(self, data: dict):
        """
        Update the configuration item with the dict provided and save it.

        :param data: dict of configuration items
        """

        self.__config.update(data)
        self.save()

    def __getattr__(self, name: str):
        """
        如果以点号方式访问属性且属性不存在于类中，尝试从配置字典中获取。
        """
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not get '{name}' in the config of the current experiment")

    def __getitem__(self, name: str):
        """
        以字典方式获取配置项的值。
        """
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not get '{name}' in the config of the current experiment")

    def __delattr__(self, name: str) -> bool:
        """
        删除配置项，如果配置项不存在,跳过

        Parameters
        ----------
        name : str
            配置项名称

        Returns
        ----------
        bool
            是否删除成功
        """
        try:
            del self.__config[name]
            return True
        except KeyError:
            return False

    def __delitem__(self, name: str) -> bool:
        """
        删除配置项，如果配置项不存在,跳过

        Parameters
        ----------
        name : str
            配置项名称

        Returns
        ----------
        bool
            是否删除成功
        """
        try:
            del self.__config[name]
            return True
        except KeyError:
            return False

    def save(self):
        """
        保存config为json，不必校验config的YAML格式，将在写入时完成校验
        """
        if not self.should_shave:
            return
        swanlog.debug("Save config to {}".format(self.__settings.get("save_path")))

        serialization_config = self.__check_config(self.__config)
        with open(self.__settings.get("save_path"), "w") as f:
            # 将config的每个key的value转换为desc和value两部分，value就是原来的value，desc是None
            # 这样做的目的是为了在web界面中显示config的内容,desc是用于描述value的
            config = {
                key: {
                    "desc": None,
                    "sort": index,
                    "value": value,
                }
                for index, (key, value) in enumerate(serialization_config.items())
            }
            yaml.dump(config, f)

    def __iter__(self):
        """
        返回配置字典的迭代器。
        """
        return iter(self.__config)

    def __len__(self):
        """
        返回配置项的数量。
        """
        return len(self.__config)

    def __str__(self):
        return str(self.__config)
