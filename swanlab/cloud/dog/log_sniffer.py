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
from ..utils import ThreadTaskABC, ThreadUtil, LogQueue
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from watchdog.utils.dirsnapshot import DirectorySnapshot, DirectorySnapshotDiff
from typing import List, Tuple, Callable
from queue import Queue
from swanlab.log import swanlog
from ..types import LogPackageType

q = Queue()


class LogSnifferHandler(FileSystemEventHandler):
    """
    日志嗅探处理器，负责处理日志信息的改动，在使用中应该对应到每一个具体类型的文件夹的监听
    由于日志信息的改动可能是十分频繁的，因此需要设置定时器，定时处理一段时间内的所有改动
    处理方式应该由传入的LogType事先定义好的处理函数决定，这里应该只负责收集日志信息
    """

    def __init__(self, watched_path: str, package_type: LogPackageType):
        """
        初始化日志嗅探处理器
        :param watched_path: 监听的路径，用作初始对照
        :param package_type: 日志包装类型，枚举类型，里面包含了不同的日志包装方式和处理方式
        """
        self.watched_path = watched_path
        self.package_type = package_type
        self.queue = LogQueue(readable=False, writable=True, new=q)
        # 如果是元数据类型，在监听之前向queue中添加上传钩子
        if package_type == LogPackageType.META:
            self.queue.put(package_type.value.before(watched_path))

    def on_modified(self, event):
        """
        在设计上所有文件都不会被删除，只会被修改，因此只需要处理文件修改事件而不需要管其他的
        """
        if event.is_directory:
            # 文件夹修改，不做处理
            # print(f"Directory modified: {event.src_path}, watched directory: {self.watched_path}")
            pass
        else:
            swanlog.debug(f"File modified: {event.src_path}, watched directory: {self.watched_path}")


class LogSnifferTask(ThreadTaskABC):
    """
    日志嗅探器，负责监听日志信息的改动，当日志信息发生改动时将改动的部分包装发送给日志聚合器
    """
    SNIFFER_TIMEOUT = 2

    def __init__(self, observer_path: List[Tuple[str, LogPackageType]]):
        """
        初始化日志嗅探器
        :param observer_path: 监听的路径和路径下对应的日志类型
        """
        self.__sniffer_queue = LogQueue(readable=True, writable=False, new=q)
        """
        日志嗅探器队列，用于存放从LogSnifferHandler中收集到的日志信息
        """
        self.__observer = Observer(timeout=self.SNIFFER_TIMEOUT)
        for watched_path, log_type in observer_path:
            self.__observer.schedule(LogSnifferHandler(watched_path, log_type),
                                     watched_path,
                                     recursive=True)
        # observer，启动！
        self.__observer.start()

    async def callback(self, u: ThreadUtil, *args):
        self.__observer.stop()

    async def task(self, u: ThreadUtil, *args):
        """
        任务执行函数，在此处收集处理的所有日志信息，解析、包装、发送给日志聚合器
        :param u: 线程工具类
        """
        u.queue.put_all(self.__sniffer_queue.get_all())
