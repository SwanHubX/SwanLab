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
from ..utils import ThreadTaskABC, ThreadUtil
from watchdog.observers import Observer
from queue import Queue
from .sniffer_queue import SnifferQueue
from .metadata_handle import MetaHandle
from typing import List
from ..task_types import UploadType
import time


class LogSnifferTask(ThreadTaskABC):
    """
    日志嗅探器，负责监听日志信息的改动，当日志信息发生改动时将改动的部分包装发送给日志聚合器
    """
    SNIFFER_TIMEOUT = 2
    __QUEUE: Queue = Queue()

    def __init__(self, meta_path: str):
        """
        初始化日志嗅探器
        :param meta_path: 元数据文件夹路径，由于目前只有元数据文件夹需要嗅探，因此只需要传入元数据文件夹路径
        后续如果有其他需要嗅探的文件夹，可以将此处改成传入handle类
        """
        self.__sniffer_queue = SnifferQueue(self.__QUEUE, readable=True, writable=False)
        """
        日志嗅探器队列，用于存放从一系列Handler中收集到的日志信息
        """
        self.__observer = Observer(timeout=self.SNIFFER_TIMEOUT)
        self.__observer.schedule(
            MetaHandle(self.__QUEUE, watched_path=meta_path),
            meta_path,
            recursive=True
        )
        # observer，启动！
        self.__observer.start()

    def callback(self, u: ThreadUtil, *args):
        # 文件事件可能会有延迟，因此需要等待一段时间
        time.sleep(self.SNIFFER_TIMEOUT)
        self.__observer.stop()
        self.pass_msg(u)

    def pass_msg(self, u: ThreadUtil):
        all_sniffer_msg: List = self.__sniffer_queue.get_all()
        if not all_sniffer_msg or len(all_sniffer_msg) == 0:
            return
        # 去重，由于现在只有files元数据文件，所以只需要针对它去重就行
        # 遍历所有的消息
        files = {UploadType.FILE: []}
        for msg in all_sniffer_msg:
            for path in msg[0]:
                if path not in files[UploadType.FILE]:
                    files[UploadType.FILE].append(path)
        new_msg = (UploadType.FILE, files[UploadType.FILE])
        u.queue.put(new_msg)

    def task(self, u: ThreadUtil, *args):
        """
        任务执行函数，在此处收集处理的所有日志信息，解析、包装、发送给日志聚合器
        :param u: 线程工具类
        """
        # 在此处完成日志信息聚合
        self.pass_msg(u)
