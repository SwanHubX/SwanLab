#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 16:41
@File: _log_sniffer.py
@IDE: pycharm
@Description:
    基于watchdog的日志资源嗅探器
"""
from .utils import LogQueue


class LogSniffer:
    """
    日志资源嗅探器，负责监控指定路径下的日志文件变化
    """

    def __init__(self, q: LogQueue):
        self.q = q
        """
        线程通信管道
        """
        self.file_path = ""
        """
        日志文件路径
        """
        self.file = None
        """
        日志文件对象
        """
        self.offset = 0
        """
        日志文件偏移量
        """
        self.buffer = ""
        """
        缓冲区
        """

    def start(self):
        """
        启动日志资源嗅探器
        """
        pass

    def stop(self):
        """
        停止日志资源嗅探器
        """
        pass

    def read(self):
        """
        读取日志文件
        """
        pass

    def parse(self):
        """
        解析日志文件
        """
        pass
