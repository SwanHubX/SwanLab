#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 19:40
@File: _log_collector.py
@IDE: pycharm
@Description:
    日志集合和上传记录器
"""

from typing import List
from .utils import ThreadUtil, ThreadTaskABC
import asyncio
from swanlab.log import swanlog
from .utils import LogQueue
from swanlab.error import UpLoadError


class LogCollectorTask(ThreadTaskABC):
    """
    日志聚合器，负责收集所有线程注册的日志信息
    并且定义日志上传接口
    """
    UPLOAD_TIME = 1
    """
    每隔多少秒上传一次日志
    """
    UPLOAD_URL = ""
    """
    日志上传的地址
    """

    def __init__(self):
        self.container: List[LogQueue.MsgType] = []
        """
        日志容器，存储从管道中获取的日志信息
        """
        self.lock = False

    async def upload(self):
        """
        将收集到的所有上传事件统一触发，上传日志信息
        所有的请求都是网络请求，因此需要异步处理，并且在此处统一
        """
        # TODO 根据日志类型进行降重
        # TODO 根据类型优先级排序，优先上传重要的日志
        tasks = [x[0].value['upload'](x[1]) for x in self.container]
        results = await asyncio.gather(*tasks)
        # 检查每一个结果
        error_tasks_index = []
        for index, result in enumerate(results):
            # 如果出现问题
            if isinstance(result, UpLoadError):
                error_tasks_index.append(index)
                continue
            elif isinstance(result, Exception):
                # 如果出现其他问题，没有办法处理，就直接跳过，但是会有警告
                swanlog.error(f"upload logs error: {result}, it might be a swanlab bug, data will be lost!")
                continue
        # 如果没有错误，清空容器
        if not len(error_tasks_index) == 0:
            return self.container.clear()
        # 如果有错误，只删除错误的部分，等待重新上传
        for index in error_tasks_index:
            self.container.pop(index)

    async def task(self, u: ThreadUtil, *args):
        """
        定时任务，定时上传日志信息
        :param u: 线程工具类
        """
        # 从管道中获取所有的日志信息，存储到self.container中
        if self.lock:
            return swanlog.debug("upload task still in progressing, passed")
        self.container.extend(u.queue.get_all())
        # print("线程" + u.name + "获取到的日志信息: ", self.container)
        if u.timer.can_run(self.UPLOAD_TIME, len(self.container) == 0):
            self.lock = True
            await self.upload()
            self.lock = False

    async def callback(self, u: ThreadUtil, *args):
        """
        回调函数，用于线程结束时的回调
        :param u: 线程工具类
        """
        print("线程" + u.name + "回调执行")
        self.container.extend(u.queue.get_all())
        return await self.upload()
