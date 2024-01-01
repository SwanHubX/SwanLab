#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    swanlab全局共用环境变量(运行时环境变量)
"""
import os
from typing import List, MutableMapping, Optional, Union

Env = Optional[MutableMapping]

# '描述' = "key"
# ---------------------------------- 基础环境变量 ----------------------------------

ROOT = "SWANLAB_ROOT"
"""命令执行目录SWANLAB_ROOT，在这个目录下寻找swanlog文件夹"""


class SwanlabConfig(object):
    """Swanlab全局配置对象"""

    def __init__(self) -> None:
        # 标志位，用于判断是否已经初始化
        self.__init = False
        # 根目录，这将决定日志输出的位置以及服务读取的位置
        self.__folder = None
        # 当前实验名称
        self.__exp_name = None
        # 当前模式，可选值: train, server; 前者代表日志记录模式，后者代表服务模式
        self.__mode = None

    def _should_initialized(func):
        """装饰器：必须在初始化完毕以后才能执行"""

        def wrapper(cls, *args, **kwargs):
            if cls.__init is False:
                raise ValueError("config has not been initialized")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def _should_added_exp(func):
        """装饰器：比如已经添加了实验"""

        def wrapper(cls, *args, **kwargs):
            if cls.__exp_name is None:
                raise ValueError("config has not add experiment")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def _should_server_mode(func):
        """装饰器：必须是server mode"""

        def wrapper(cls, *args, **kwargs):
            if cls.__mode != "server":
                raise ValueError(f"{func.__name__} is only available in server mode")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def _should_train_mode(func):
        """装饰器：必须是train mode"""

        def wrapper(cls, *args, **kwargs):
            if cls.__mode != "train":
                raise ValueError(f"{func.__name__} is only available in train mode")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def init(self, root: str, mode: str):
        """初始化配置对象"""
        if self.__init:
            # FIXME 可能需要输出一下，说当前已经初始化过了
            return
        self.__folder = root
        if mode not in ["train", "server"]:
            raise ValueError("mode must be train or server")
        self.__mode = mode
        self.__init = True

    def add_exp(self, exp_name: str):
        if self.__exp_name is not None and self.__mode == "train":
            raise ValueError("config has been added experiment in train mode")
        self.__exp_name = exp_name

    @property
    @_should_initialized
    def isTrain(self) -> str:
        """当前模式是否为训练模式"""
        return self.__mode == "train"

    @staticmethod
    def getcwd() -> str:
        """当前程序运行路径，不包括文件名"""
        return os.getcwd()

    @property
    @_should_initialized
    @_should_added_exp
    def exp_name(self):
        return self.__exp_name

    @property
    @_should_initialized
    def root(self) -> str:
        """项目输出根目录，必须先被初始化"""
        r = os.path.join(self.__folder, "swanlog")
        if not os.path.exists(r):
            os.mkdir(r)
        return r

    @property
    @_should_initialized
    def project(self) -> str:
        """项目配置文件路径，必须是训练模式"""
        return os.path.join(self.root, "project.json")

    @property
    @_should_initialized
    def output(self) -> str:
        """服务日志输出文件路径或者训练时swanlab的日志输出文件路径"""
        if self.__mode == "train":
            return os.path.join(self.exp_folder, "output.log")
        else:
            #
            return os.path.join(self.root, "output.log")

    @property
    @_should_initialized
    @_should_train_mode
    @_should_added_exp
    def exp_folder(self) -> str:
        """实验存储路径"""
        return os.path.join(self.root, self.__exp_name)

    @property
    @_should_initialized
    @_should_train_mode
    @_should_added_exp
    def logs_folder(self) -> str:
        """日志输出根目录，必须是训练模式"""
        return os.path.join(self.root, self.__exp_name, "logs")

    @property
    @_should_initialized
    @_should_train_mode
    @_should_added_exp
    def chart(self) -> str:
        """表格路径"""
        return os.path.join(self.root, self.__exp_name, "chart.json")

    @property
    @_should_initialized
    @_should_train_mode
    @_should_added_exp
    def console_folder(self) -> str:
        """终端监听文件根目录，必须是训练模式"""
        return os.path.join(self.root, self.__exp_name, "console")

    @property
    @_should_initialized
    @_should_train_mode
    @_should_added_exp
    def error(self) -> str:
        """终端错误日志打印路径"""
        return os.path.join(self.root, self.__exp_name, "console", "error.log")


swc = SwanlabConfig()


# ---------------------------------- 定义变量访问方法 ----------------------------------


def get_runtime_root(env: Optional[Env] = None) -> Optional[str]:
    """获取运行时根路径

    Parameters
    ----------
    env : Optional[Env], optional
        环境变量map, by default None

    Returns
    -------
    Optional[str]
        根路径
    """
    default: Optional[str] = os.getcwd()
    if env is None:
        env = os.environ
    # ROOT拿到的应该是一个绝对路径
    return env.get(ROOT, default=default)


def get_swanlog_dir() -> Optional[str]:
    """获取swanlog路径，这是一个计算变量，
    通过`get_runtime_root()`返回值得到

    Returns
    -------
    Optional[str]
        swanlog目录路径
    """
    return os.path.join(get_runtime_root(), "swanlog")


# TODO 后续改为数据库路径
def get_runtime_project() -> Optional[str]:
    """获取运行时项目配置，这是一个计算变量，
    通过`get_swanlog_dir()`返回值得到

    Returns
    -------
    Optional[str]
        项目配置文件路径
    """
    return os.path.join(get_swanlog_dir(), "project.json")
