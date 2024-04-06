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
from .files_types import FileType


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

    @staticmethod
    def report_known_error(errors: List[UpLoadError]):
        """
        上报错误信息
        :param errors: 错误信息列表
        """
        if len(errors) == 0:
            return
        # 去重
        errors = list(set(errors))
        for error in errors:
            swanlog.__getattribute__(error.log_level)(error.message)

    async def upload(self):
        """
        将收集到的所有上传事件统一触发，上传日志信息
        所有的请求都是网络请求，因此需要异步处理，并且在此处统一
        """
        # 根据日志类型进行降重
        tasks_dict = {x: [] for x in FileType}
        for msg in self.container:
            tasks_dict[msg[0]].extend(msg[1])
        # 此时应该只剩下最多 FileType内部枚举个数 个任务
        # 媒体类型需要首先上传，后面几个可以一起上传
        # 检查每一个上传结果
        success_tasks_type = []
        # 已知错误列表
        known_errors = []
        if len(tasks_dict[FileType.MEDIA]):
            media_result = await FileType.MEDIA.value['upload'](tasks_dict[FileType.MEDIA])
            if isinstance(media_result, UpLoadError):
                known_errors.append(media_result)
            elif isinstance(media_result, Exception):
                swanlog.error(f"upload logs error: {media_result}, it might be a swanlab bug, data will be lost!")
            success_tasks_type.append(FileType.MEDIA)
        # 上传任务
        tasks_key_list = [key for key in tasks_dict if len(tasks_dict[key]) > 0 and key != FileType.MEDIA]
        tasks = [x.value['upload'](tasks_dict[x]) for x in tasks_key_list]
        results = await asyncio.gather(*tasks)
        for index, result in enumerate(results):
            # 如果出现已知问题
            if isinstance(result, UpLoadError):
                known_errors.append(result)
                continue
            # 如果出现其他问题，没有办法处理，就直接跳过，但是会有警告
            elif isinstance(result, Exception):
                swanlog.error(f"upload logs error: {result}, it might be a swanlab bug, data will be lost!")
                continue
            # 标记所有已经成功的任务
            success_tasks_type.append(tasks_key_list[index])
        # 将不成功的任务加入原本的容器中
        self.container = [(x, tasks_dict[x]) for x in tasks_dict if x not in success_tasks_type]
        self.report_known_error(known_errors)

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
