#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/4 14:11
@File: log_sniffer.py
@IDE: pycharm
@Description:
    日志嗅探器
"""
from enum import Enum
from ..utils import ThreadTaskABC, ThreadUtil
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


# from swanlab.db


class LogPackageSeed:
    """
    日志包装种子，列举不同的日志包装方式
    """
    pass


class LogType(Enum):
    """
    日志类型枚举，涉及不同的日志包装方式
    """
    pass


class LogSnifferHandler(FileSystemEventHandler):
    def on_modified(self, event):
        pass


class LogSnifferTask(ThreadTaskABC):
    """
    日志嗅探器，负责监听日志信息的改动，当日志信息发生改动时将改动的部分包装发送给日志聚合器
    """

    async def callback(self, u: ThreadUtil, *args):
        pass

    async def task(self, u: ThreadUtil, **kwargs):
        pass

    def __init__(self):
        pass

    def run(self):
        pass
