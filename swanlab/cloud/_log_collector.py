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
from .utils import ThreadUtil
import asyncio


class LogCollector:
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
        self.container: List = []
        """
        日志容器，存储从管道中获取的日志信息
        """

    async def upload(self):
        """
        上传日志信息，异步上传
        """
        # TODO 上传日志信息
        # 现在暂时在运行时在本地文件中写入追加的日志信息
        with open("log.txt", "a") as f:
            print("开始上传日志信息:", self.container)
            for msg in self.container:
                f.write(msg + "\n")
            f.write('\n\n')
        # 假设上传时间为1秒
        await asyncio.sleep(1)

    async def task(self, u: ThreadUtil):
        """
        定时任务，定时上传日志信息
        :param u: 线程工具类
        """
        # 从管道中获取所有的日志信息，存储到self.container中
        self.container.extend(u.queue.get_all())
        print("线程" + u.name + "获取到的日志信息: ", self.container)
        if u.timer.can_run(self.UPLOAD_TIME, len(self.container) == 0):
            await self.upload()
            # 清除容器内容
            self.container.clear()

    async def callback(self, u: ThreadUtil):
        """
        回调函数，用于线程结束时的回调
        :param u: 线程工具类
        """
        print("线程" + u.name + "回调执行")
        self.container.extend(u.queue.get_all())
        return await self.upload()
