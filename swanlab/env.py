#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    环境变量配置，用于配置一些全局变量，如存储路径等内容
"""
import os
import mimetypes
from functools import wraps

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
因此，需要指定文件路径与文件名，用于后端服务的响应（这在下面的路由配置中也有说明）
"""
# 注册静态文件路径
mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")
# 静态文件路径
FILEPATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")


class SwanlabConfig(object):
    """Swanlab全局配置对象"""

    # 标志位，用于判断是否已经初始化
    __init = False

    def __init__(self) -> None:
        # 根目录，这将决定日志输出的位置以及服务读取的位置
        self.__folder = None
        # 当前实验名称
        self.__exp_name = None
        # 当前模式，可选值: train, server; 前者代表日志记录模式，后者代表服务模式
        self.__mode = None

    def __should_initialized(func):
        """装饰器：必须在初始化完毕以后才能执行"""
        wraps(func)

        def wrapper(cls, *args, **kwargs):
            if cls.__init is False:
                raise ValueError("config has not been initialized")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def __should_added_exp(func):
        """装饰器：比如已经添加了实验"""
        wraps(func)

        def wrapper(cls, *args, **kwargs):
            if cls.__exp_name is None:
                raise ValueError("config has not add experiment")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def __should_server_mode(func):
        """装饰器：必须是server mode"""
        wraps(func)

        def wrapper(cls, *args, **kwargs):
            if cls.__mode != "server":
                raise ValueError(f"{func.__name__} is only available in server mode")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def __should_train_mode(func):
        """装饰器：必须是train mode"""
        wraps(func)

        def wrapper(cls, *args, **kwargs):
            if cls.__mode != "train":
                raise ValueError(f"{func.__name__} is only available in train mode")
            result = func(cls, *args, **kwargs)
            return result

        return wrapper

    def init(self, root: str, mode: str):
        """初始化配置对象"""
        if self.__init:
            raise ValueError("config has been initialized")
        self.__folder = root
        if mode not in ["train", "server"]:
            raise ValueError("mode must be train or server")
        self.__mode = mode
        self.__init = True

    def add_exp(self, exp_name: str):
        self.__exp_name = exp_name

    @staticmethod
    def getcwd() -> str:
        """当前程序运行路径，不包括文件名"""
        return os.getcwd()

    @property
    @__should_initialized
    def root(self) -> str:
        """项目输出根目录，必须先被初始化"""
        return os.path.join(self.__folder, "swanlog")

    @property
    @__should_initialized
    def project(self) -> str:
        """项目配置文件路径，必须是训练模式"""
        return os.path.join(self.root, "project.json")

    @property
    @__should_initialized
    @__should_server_mode
    def output(self) -> str:
        """服务日志输出文件夹"""
        return os.path.join(self.root, "output.log")

    @property
    @__should_initialized
    @__should_added_exp
    def exp_folder(self) -> str:
        """实验存储路径"""
        return os.path.join(self.root, self.__exp_name)

    @property
    @__should_initialized
    @__should_train_mode
    @__should_added_exp
    def logs_folder(self) -> str:
        """日志输出根目录，必须是训练模式"""
        return os.path.join(self.root, self.__exp_name, "logs")

    @property
    @__should_initialized
    @__should_train_mode
    @__should_added_exp
    def chart(self) -> str:
        """表格路径"""
        return os.path.join(self.root, self.__exp_name, "chart.json")

    @property
    @__should_initialized
    @__should_train_mode
    @__should_added_exp
    def console_folder(self) -> str:
        """终端监听文件根目录，必须是训练模式"""
        return os.path.join(self.root, self.__exp_name, "console")

    @property
    @__should_initialized
    @__should_train_mode
    @__should_added_exp
    def output_folder(self) -> str:
        """项目日志输出位置文件夹"""
        return os.path.join(self.root, self.__exp_name, "output")


swc = SwanlabConfig()


if __name__ == "__main__":
    swc.init("test", "server")
    print(swc.root)
    print(swc.output)
    print(swc.console_folder)
