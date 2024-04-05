#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/4 14:11
@File: log_sniffer.py
@IDE: pycharm
@Description:
    日志嗅探器
    嗅探不做上传操作，只做采集操作，将采集到的日志、差异信息发送给日志聚合器
"""
from ..utils import ThreadTaskABC, ThreadUtil, LogQueueABC
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from watchdog.utils.dirsnapshot import DirectorySnapshot, DirectorySnapshotDiff
from typing import List, Tuple, Dict
from queue import Queue
import asyncio

q = Queue()


class SnifferQueue(LogQueueABC):

    def __init__(self, readable: bool = True, writable: bool = True):
        super().__init__(q, readable, writable)


class LogSnifferHandler(FileSystemEventHandler):
    """
    日志嗅探处理器，负责处理日志信息的改动，在使用中应该对应到每一个具体类型的文件夹的监听
    由于日志信息的改动可能是十分频繁的，因此需要设置定时器，定时处理一段时间内的所有改动
    处理方式应该由传入的LogType事先定义好的处理函数决定，这里应该只负责收集日志信息
    Watchdog 默认情况下使用两个线程。一个线程用于监视文件系统事件，另一个线程用于处理这些事件并触发相应的回调函数
    因此在这里不需要考虑多线程问题
    """

    def __init__(self, watched_path: str):
        """
        初始化日志嗅探处理器
        :param watched_path: 监听的路径，用作初始对照
        """
        self.watched_path = watched_path
        self.queue = SnifferQueue(readable=False, writable=True)

    def on_modified(self, event):
        """
        在设计上所有文件都不会被删除，只会被修改，因此只需要处理文件修改事件而不需要管其他的
        """
        if event.is_directory:
            # 文件夹修改，不做处理
            # print(f"Directory modified: {event.src_path}, watched directory: {self.watched_path}")
            pass
        else:
            print(f"File modified: {event.src_path}, watched directory: {self.watched_path}")


class LogSnifferTask(ThreadTaskABC):
    """
    日志嗅探器，负责监听日志信息的改动，当日志信息发生改动时将改动的部分包装发送给日志聚合器
    """
    SNIFFER_TIMEOUT = 2

    def __init__(self, observer_paths: List[Tuple[str,]]):
        """
        初始化日志嗅探器
        :param observer_paths: 监听的路径和路径下对应的日志类型
        """
        self.__sniffer_queue = SnifferQueue(readable=True, writable=False)
        """
        日志嗅探器队列，用于存放从LogSnifferHandler中收集到的日志信息
        """
        self.__observer = Observer(timeout=self.SNIFFER_TIMEOUT)
        for watched_path, in observer_paths:
            self.__observer.schedule(LogSnifferHandler(watched_path),
                                     watched_path,
                                     recursive=True)
        # observer，启动！
        self.__observer.start()

    async def callback(self, u: ThreadUtil, *args):
        # 文件事件可能会有延迟，因此需要等待一段时间
        # await asyncio.sleep(1.5)
        self.__observer.stop()

    async def task(self, u: ThreadUtil, *args):
        """
        任务执行函数，在此处收集处理的所有日志信息，解析、包装、发送给日志聚合器
        :param u: 线程工具类
        """
        # 在此处完成日志信息聚合
        print("日志嗅探器开始执行")
