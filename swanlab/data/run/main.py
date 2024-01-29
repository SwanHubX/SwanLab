#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:58:58
@File: swanlab/data/run/main.py
@IDE: vscode
@Description:
    在此处定义SwanLabRun类并导出
"""
from typing import Any
from ..settings import SwanDataSettings
from ...log import register, swanlog
from ..system import get_system_info, get_requirements
from .utils import (
    check_exp_name_format,
    check_desc_format,
    get_a_lock,
    json_serializable,
)
from datetime import datetime
import os, time
import random
import ujson
from .exp import SwanLabExp
from collections.abc import Mapping
from .db import Experiment, ExistedError
from typing import Tuple
import yaml
import argparse


def need_inited(func):
    """装饰器，用于检查是否已经初始化"""

    def wrapper(self, *args, **kwargs):
        if not self._inited:
            raise RuntimeError("You must call swanlab.init() before using swanlab.log")
        return func(self, *args, **kwargs)

    return wrapper


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

    def __init__(self, config: dict, settings: SwanDataSettings = None):
        """
        实例化配置类，如果settings不为None，说明是通过swanlab.init调用的，否则是通过swanlab.config调用的

        Parameters
        ----------
        settings : SwanDataSettings, optional
            运行时设置
        """
        self.__config.update(self.__check_config(config))
        self.__settings["save_path"] = settings.config_path if settings is not None else None
        if self._inited:
            self.__save()

    def __check_config(self, config: dict) -> dict:
        """
        检查配置是否合法，确保它可以被 JSON/YAML 序列化。
        如果传入的是 argparse.Namespace 类型，会先转换为字典。
        """
        if config is None:
            return {}
        # config必须可以被json序列化
        try:
            if isinstance(config, argparse.Namespace):
                config = vars(config)
            # 将config转换为json序列化的dict
            config = json_serializable(dict(config))
            # 尝试序列化，如果还是失败就退出
            yaml.dump(config)
        except:
            raise TypeError(f"config: {config} is not a valid dict, which can be json serialized")
        return config

    def __check_private(self, name: str):
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

    @need_inited
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
        self.__check_private(name)
        # 设置属性，并判断是否已经初始化，如果是，则调用保存方法
        self.__dict__[name] = value
        # 同步到配置字典
        self.__config[name] = value
        self.__save()

    @need_inited
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
        self.__check_private(name)
        self.__config[name] = value
        self.__save()

    @need_inited
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
        self.__check_private(name)
        self.__config[name] = value
        self.__save()

    @need_inited
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
            self.__save()
            return True
        except KeyError:
            return False

    @need_inited
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

    @need_inited
    def __getattr__(self, name: str):
        """
        如果以点号方式访问属性且属性不存在于类中，尝试从配置字典中获取。
        """
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not get '{name}' in the config of the current experiment")

    @need_inited
    def __getitem__(self, name: str):
        """
        以字典方式获取配置项的值。
        """
        try:
            return self.__config[name]
        except KeyError:
            raise AttributeError(f"You have not get '{name}' in the config of the current experiment")

    @need_inited
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

    @need_inited
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

    def __save(self):
        """
        保存config为json，不必校验config的YAML格式，将在写入时完成校验
        """
        swanlog.debug("Save config to {}".format(self.__settings.get("save_path")))
        with get_a_lock(self.__settings.get("save_path"), "w") as f:
            # 将config的每个key的value转换为desc和value两部分，value就是原来的value，desc是None
            # 这样做的目的是为了在web界面中显示config的内容,desc是用于描述value的
            config = {
                key: {
                    "desc": None,
                    "sort": index,
                    "value": value,
                }
                for index, (key, value) in enumerate(self.__config.items())
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


class SwanLabRun:
    """
    The SwanLabRun class is used for logging during a single experiment.
    There should be only one instance of the SwanLabRun class for each experiment.
    """

    def __init__(
        self,
        experiment_name: str = None,
        description: str = None,
        config: dict = None,
        log_level: str = None,
        suffix: str = None,
        loggings: bool = False,
    ):
        """
        Initializing the SwanLabRun class involves configuring the settings and initiating other logging processes.

        Parameters
        ----------
        experiment_name : str, optional
            实验名称，实验名称应该唯一，由0-9，a-z，A-Z，" ","_","-","/"组成
            如果不提供此参数(为None)，SwanLab将自动生成一个实验名称
        description : str, optional
            实验描述，用于对当前实验进行更详细的介绍或标注
            如果不提供此参数(为None)，可以在web界面中进行修改,这意味着必须在此改为空字符串""
        config : dict, optional
            实验参数配置，可以在web界面中显示，如学习率、batch size等
            不需要做任何限制，但必须是字典类型，可被json序列化，否则会报错
        log_level : str, optional
            当前实验的日志等级，默认为'info'，可以从'debug'、'info'、'warning'、'error'、'critical'中选择
            不区分大小写，如果不提供此参数(为None)，则默认为'info'
            如果提供的日志等级不在上述范围内，默认改为info
        suffix : str, optional
            实验名称后缀，用于区分同名实验，格式为yyyy-mm-dd_HH-MM-SS
            如果不提供此参数(为None)，不会添加后缀
        """
        # ---------------------------------- 初始化类内参数 ----------------------------------
        # 生成一个唯一的id，随机生成一个8位的16进制字符串，小写
        id = hex(random.randint(0, 2**32 - 1))[2:].zfill(8)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.__run_id = "run-{}-{}".format(timestamp, id)
        # 初始化配置
        self.__settings = SwanDataSettings(run_id=self.__run_id)
        # ---------------------------------- 初始化日志记录器 ----------------------------------
        # output、console_dir等内容不依赖于实验名称的设置
        register(self.__settings.output_path, self.__settings.console_dir, enable_logging=loggings)
        # 初始化日志等级
        level = self.__check_log_level(log_level)
        swanlog.setLevel(level)

        # ---------------------------------- 初始化配置 ----------------------------------
        # 给外部1个config
        self.__config = SwanLabConfig(config, self.__settings)
        # ---------------------------------- 注册实验 ----------------------------------
        # 校验描述格式
        description = self.__check_description(description)
        self.__exp: SwanLabExp = self.__register_exp(
            experiment_name,
            description,
            suffix,
        )
        # 实验状态标记，如果status不为0，则无法再次调用log方法
        self.__status = 0

    @property
    def settings(self) -> SwanDataSettings:
        """
        This property allows you to access the 'settings' content passed through `init`,
        and runtime settings can not be modified.
        """
        return self.__settings

    @property
    def config(self):
        """
        This property allows you to access the 'config' content passed through `init`,
        and allows you to modify it. The latest configuration after each modification
        will be synchronized to the corresponding path by Swanlab. If you have
        enabled the web service, you will notice the changes after refreshing.
        """
        return self.__config

    def log(self, data: dict, step: int = None):
        """
        Log a row of data to the current run.
        Unlike `swanlab.log`, this api will be called directly based on the SwanRun instance, removing the initialization process.
        Of course, after you call the success/fail method, this method will be banned from calling.

        Parameters
        ----------
        data : Dict[str, DataType]
            Data must be a dict.
            The key must be a string with 0-9, a-z, A-Z, " ", "_", "-", "/".
            The value must be a `float`, `float convertible object`, `int` or `swanlab.data.BaseType`.
        step : int, optional
            The step number of the current data, if not provided, it will be automatically incremented.
            If step is duplicated, the data will be ignored.
        """
        if self.__status != 0:
            raise RuntimeError("After experiment finished, you can no longer log data to the current experiment")

        if not isinstance(data, dict):
            return swanlog.error(
                "log data must be a dict, but got {}, SwanLab will ignore records it.".format(type(data))
            )
        # 检查step的类型
        if step is not None and (not isinstance(step, int) or step < 0):
            swanlog.error(
                "'step' must be an integer not less than zero, but got {}, SwanLab will automatically set step".format(
                    step
                )
            )
            step = None
        # 遍历data，记录data
        for key in data:
            # 遍历字典的key，记录到本地文件中
            d = data[key]
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            self.__exp.add(key=key, data=d, step=step)

    def success(self):
        """标记实验成功"""
        self.__set_exp_status(1)

    def fail(self):
        """标记实验失败"""
        self.__set_exp_status(-1)

    def __str__(self) -> str:
        """此类的字符串表示"""
        return self.__run_id

    def __set_exp_status(self, status: int):
        """更新实验状态

        Parameters
        ----------
        status : int
            实验状态，1代表成功，-1代表失败，0代表未完成，其他值会被自动设置为-1
        """
        if status not in [1, -1, 0]:
            swanlog.warning("Invalid status when set, status must be 1, -1 or 0, but got {}".format(status))
            swanlog.warning("SwanLab will set status to -1")
            status = -1
        self.__status = status
        self.__exp.db.status = status
        self.__exp.db.save()

    def __get_exp_name(self, experiment_name: str = None, suffix: str = None) -> Tuple[str, str]:
        """
        预处理实验名称，如果实验名称过长，截断

        Parameters
        ----------
        experiment_name : str
            实验名称
        suffix : str
            实验名称后缀添加方式，可以为None、"timestamp"，前者代表不添加后缀，后者代表添加时间戳后缀

        Returns
        ----------
        experiment_name : str
            校验后的实验名称
        exp_name : str
            最终的实验名称
        """
        # ---------------------------------- 校验实验名称 ----------------------------------
        experiment_name = "exp" if experiment_name is None else experiment_name
        # 校验实验名称
        experiment_name_checked = check_exp_name_format(experiment_name)
        # 如果前后长度不一样，说明实验名称被截断了，提醒
        if len(experiment_name_checked) != len(experiment_name):
            swanlog.warning("The experiment name you provided is not valid, it has been truncated automatically.")
        # 为实验名称添加后缀，格式为yyyy-mm-dd_HH-MM-SS
        if suffix is None:
            return experiment_name_checked, experiment_name
        # 如果suffix不是timestamp，提醒，自动改为timestamp
        if suffix.lower() != "timestamp":
            suffix = "timestamp"
            swanlog.warning(f"The suffix you provided is not valid, it has been set to {suffix}.")
        # ---------------------------------- 自动添加后缀 ----------------------------------
        # 现在只有一种后缀，即timestamp，所以直接添加就行
        exp_name = "{}_{}".format(experiment_name_checked, datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
        return experiment_name_checked, exp_name

    def __register_exp(
        self,
        experiment_name: str,
        description: str = None,
        suffix: str = None,
    ) -> SwanLabExp:
        """
        注册实验，将实验配置写入数据库中，完成实验配置的初始化
        """
        # 这个循环的目的是如果创建失败则等零点五秒重新生成后缀重新创建，直到创建成功
        # 但是由于需要考虑suffix为none不生成后缀的情况，所以需要在except中判断一下
        while True:
            experiment_name, exp_name = self.__get_exp_name(experiment_name, suffix)
            try:
                # 获得数据库实例
                exp = Experiment.create(name=exp_name, run_id=self.__run_id, description=description)
                break
            except ExistedError:
                if suffix is None:
                    raise ExistedError(f"Experiment {exp_name} has existed, please try another name.")
                # 如果suffix不为None，说明是自动生成的后缀，需要重新生成后缀
                swanlog.debug(f"Experiment {exp_name} has existed, try another name...")
                time.sleep(0.5)
                continue
        self.__settings.exp_name = exp_name
        # 实验创建成功，执行一些记录操作
        self.__record_exp_config()  # 记录实验配置
        # 打印信息

        swanlog.info(f"Experiment {experiment_name} has been registered.")
        return SwanLabExp(self.__settings, exp.id, exp=exp)

    def __check_log_level(self, log_level: str) -> str:
        """检查日志等级是否合法"""
        valids = ["debug", "info", "warning", "error", "critical"]
        if log_level is None:
            return "info"
        elif log_level.lower() in valids:
            return log_level.lower()
        else:
            swanlog.warning(f"The log level you provided is not valid, it has been set to {log_level}.")
            return "info"

    def __check_description(self, description: str) -> str:
        """检查实验描述是否合法"""
        if description is None:
            return ""
        desc = check_desc_format(description)
        if desc != description:
            swanlog.warning("The description has been truncated automatically.")
        return desc

    def __record_exp_config(self):
        """创建实验配置目录 files
        - 创建 files 目录
        - 将实验环境写入 files/swanlab-metadata.json 中
        - 将实验依赖写入 files/requirements.txt 中
        """
        requirements_path = self.__settings.requirements_path
        metadata_path = self.__settings.metadata_path
        # 将实验依赖存入 requirements.txt
        with open(requirements_path, "w") as f:
            f.write(get_requirements())
        # 将实验环境(硬件信息、git信息等等)存入 swanlab-metadata.json
        with open(metadata_path, "w") as f:
            ujson.dump(get_system_info(), f)
