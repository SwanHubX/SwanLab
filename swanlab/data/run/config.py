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
from swanlab.log import swanlog
import datetime
import math
from typing import Callable, Optional
from .callback import RuntimeInfo
from swanlab.data.modules import Line
import re
import json


def json_serializable(obj: dict):
    """
    将传入的字典转换为JSON可序列化格式。
    """
    # 如果对象是基本类型，则直接返回
    if isinstance(obj, (int, float, str, bool, type(None))):
        if isinstance(obj, float) and math.isnan(obj):
            return Line.nan
        if isinstance(obj, float) and math.isinf(obj):
            return Line.inf
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


def third_party_config_process(data) -> dict:
    """
    对于一些特殊的第三方库的处理，例如omegaconf
    """
    # 如果是omegaconf的DictConfig，则转换为字典
    try:
        import omegaconf  # noqa
        if isinstance(data, omegaconf.DictConfig):
            return omegaconf.OmegaConf.to_container(data, resolve=True, throw_on_missing=True)
        else:
            raise TypeError
    except ImportError:
        pass

    # 如果是argparse的Namespace，则转换为字典
    if isinstance(data, argparse.Namespace):
        return vars(data)
    else:
        raise TypeError


def parse(config: Mapping) -> dict:
    """
    Check the configuration item and convert it to a JSON serializable format.
    """
    if config is None:
        return {}
    # 1. 第三方配置类型判断与转换
    try:
        return third_party_config_process(config)
    except ImportError:
        pass
    except TypeError:
        pass
    # 2. 将config转换为可被json序列化的字典
    try:
        return json_serializable(dict(config))
    except Exception:  # noqa
        pass
    # 3. 尝试序列化，序列化成功直接返回
    try:
        return json.loads(json.dumps(config))
    except Exception as e:  # noqa
        # 还失败就没办法了，👋
        raise TypeError(f"config: {config} is not a json serialized dict, error: {e}")


class SwanLabConfig(Mapping):
    """
    The SwanConfig class is used for realize the invocation method of `run.config.lr`.

    Attention:
    The configuration item must be JSON serializable; Cannot set private attributes by `.__xxx`.
    """
    __dict__ = {
        "_SwanLabConfig__config": {},
        "_SwanLabConfig__on_setter": None,
    }

    def __init__(self, config: Mapping = None, on_setter: Optional[Callable[[RuntimeInfo], Any]] = None):
        """
        实例化配置类，如果settings不为None，说明是通过swanlab.init调用的，否则是通过swanlab.config调用的
        """
        if config is not None:
            self.__config.update(parse(config))
        # 每一个实例有自己的config
        self.__config = {}
        self.__on_setter = on_setter

    @staticmethod
    def __fmt_config(config: dict):
        """
        格式化config，值改为value字段，增加desc和sort字段
        """
        # 遍历每一个配置项，值改为value
        sort = 0
        for key, value in config.items():
            config[key] = {"value": value, "desc": "", "sort": sort}
            sort += 1

    def __save(self):
        """
        保存config为dict
        """
        if not self.__on_setter:
            return
        try:
            # 深度拷贝一次，防止引用传递
            data = yaml.load(yaml.dump(self.__config), Loader=yaml.FullLoader)
        except Exception as e:
            swanlog.error(f"Error occurred when saving config: {e}")
            return
        # 遍历每一个配置项，值改为value，如果是字典，则递归调用
        self.__fmt_config(data)
        r = RuntimeInfo(config=self.__config)
        self.__on_setter(r)

    # ---------------------------------- 实现对象风格 ----------------------------------

    def __delattr__(self, name: str):
        """
        删除配置项，如果配置项不存在
        """
        # _*__正则开头的属性不允许删除
        if re.match(r"_.*__", name):
            raise AttributeError(f"Attribute '{name}' is private and cannot be deleted")
        try:
            del self.__config[name]
            self.__save()
        except KeyError:
            raise AttributeError(f"You have not deleted '{name}' in the config of the current experiment")

    def __getattr__(self, name: str):
        """
        如果以点号方式访问属性且属性不存在于类中，尝试从配置字典中获取。
        """
        # 如果self.__dict__中有name属性，则返回
        try:
            return self.__dict__[name]
        except KeyError:
            pass
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not get '{name}' in the config of the current experiment")

    def __setattr__(self, name: str, value: Any) -> None:
        """
        Custom setter attribute, user can not set private attributes.
        """
        name = str(name)
        if name in self.__dict__:
            self.__dict__[name] = value
            return
        # _*__正则开头的属性不允许设置
        if re.match(r"_.*__", name):
            raise AttributeError(f"Attribute '{name}' is private and cannot be set")
        # 否则应该设置到配置字典中
        self.__config[name] = parse(value)
        self.__save()

    # ---------------------------------- 实现字典风格 ----------------------------------

    def get(self, name: str, default=None):
        """
        Get the value of a configuration item. If the item does not exist, raise AttributeError.
        """
        try:
            return self.__config[name]
        except KeyError:
            return default

    def __delitem__(self, name: str):
        """
        删除配置项，如果配置项不存在,跳过
        """
        try:
            del self.__config[name]
            self.__save()
        except KeyError:
            raise KeyError(f"You have not set '{name}' in the config of the current experiment when deleting")

    def __getitem__(self, name: str):
        """
        以字典方式获取配置项的值
        """
        # 如果self.__dict__中有name属性，则返回
        # 以_SwanLabConfig__开头，删除
        if name.startswith("_SwanLabConfig__"):
            name = name[15:]
        try:
            return self.__config[name]
        except KeyError:
            raise KeyError(f"You have not get '{name}' in the config of the current experiment")

    def __setitem__(self, name: str, value: Any) -> None:
        """
        Set the value of a configuration item. If the item does not exist, create it.
        User are not allowed to set private attributes.
        """
        name = str(name)
        self.__config[name] = parse(value)
        self.__save()

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

    # ---------------------------------- 其他函数 ----------------------------------

    def set(self, name: str, value: Any):
        """
        Explicitly set the value of a configuration item and save it.
        Private attributes are not allowed to be set.
        """
        name = str(name)
        self.__config[name] = parse(value)
        self.__save()

    def pop(self, name: str):
        """
        Delete a configuration item; if the item does not exist, skip.
        """
        try:
            t = self.__config[name]
            del self.__config[name]
            self.__save()
            return t
        except KeyError:
            return None

    def update(self, data: dict):
        """
        Update the configuration item with the dict provided and save it.
        :param data: dict of configuration items
        """
        self.__config.update(parse(data))
        self.__save()

    def clean(self):
        """
        Clean the configuration.
        Attention: This method will reset the instance and instance will not automatically save the configuration.
        """
        self.__config.clear()
        self.__on_setter = None
