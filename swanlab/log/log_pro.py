"""
@author: cunyue
@file: log_pro.py
@time: 2025/5/15 18:35
@description: 标准输出、标准错误流拦截代理，支持外界设置/取消回调，基础作用为输出日志
"""

import sys
from typing import Callable

from swankit.log import SwanLabSharedLog

from swanlab.swanlab_settings import get_settings

SwanLogCallback = Callable[]

class SwanLog(SwanLabSharedLog):
    """
    swanlab 日志类
    继承自 SwanLabSharedLog 的同时增加标准输出、标准错误留拦截代理功能
    """

    def __init__(self, name=__name__.lower(), level="info"):
        super().__init__(name=name, level=level)
        self.__log_state_manager = LogStateManager()


class LogStateManager:
    """
    日志状态管理器，统一保存一些日志状态
    """

    def __init__(self):
        # 保存原始的标准输出和标准错误留
        self.__origin_stdout_write = None
        self.__origin_stderr_write = None
        # 当前已经代理的输出行数
        self.__epoch = 0
        # 代理缓冲区
        self.__buffer = ""
        # 上传到云端的最大长度
        self.__max_upload_len = get_settings().max_log_length

    @property
    def __proxied(self):
        """
        判断是否已经开启代理了
        :return: bool
        """
        return self.__origin_stderr_write is not None and self.__origin_stdout_write is not None

    def create_write_handler(self, callback: Callable, is_error: bool = False) -> Callable:
        """
        生成一个新的写入处理器，可选写入输出流还是错误流
        :return: Callable
        """
        origin_write_handler = self.__origin_stderr_write if is_error else self.__origin_stdout_write

        # 生成新的写入处理器
        def write_handler(message: str):
            # 1. 原模原样写入原本的输出流中
            try:
                origin_write_handler(message)
            except UnicodeEncodeError:
                # 遇到编码问题，直接pass，此时表现为终端不输出
                pass
            except ValueError as e:
                # 遇到文件已关闭问题，直接pass，此时表现为终端不输出
                if "I/O operation on closed file" in str(e):
                    pass
            # 2. 自己代理一些状态
            self.__epoch += 1

        return write_handler

    def start_proxy(self, callback: Callable):
        """
        启动代理
        :param callback: 回调函数列表
        :return: None
        """
        if self.__proxied:
            raise RuntimeError("Std Proxy is already started")
        self.__origin_stdout_write = sys.stdout.write
        self.__origin_stderr_write = sys.stderr.write
        self.__epoch = 0
        self.__buffer = ""
        sys.stdout.write = self.create_write_handler(callback, is_error=False)
        sys.stderr.write = self.create_write_handler(callback, is_error=True)

    def stop_proxy(self):
        """
        停止代理
        :return: None
        """
        if not self.__proxied:
            raise RuntimeError("Std Proxy is not started")
        sys.stdout.write = self.__origin_stdout_write
        sys.stderr.write = self.__origin_stderr_write
        self.__origin_stdout_write = None
        self.__origin_stderr_write = None
        self.__epoch = 0
        self.__buffer = ""
