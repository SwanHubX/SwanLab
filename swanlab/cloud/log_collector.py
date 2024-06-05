#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 19:40
@File: log_collector.py
@IDE: pycharm
@Description:
    日志集合和上传记录器
"""
import time
from typing import List
from .utils import ThreadUtil, ThreadTaskABC
from swanlab.log import swanlog
from .utils import LogQueue
from swanlab.error import SyncError
from .task_types import UploadType


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

    def __init__(self, upload_type=UploadType):
        self.container: List[LogQueue.MsgType] = []
        """
        日志容器，存储从管道中获取的日志信息
        """
        self.lock = False
        self.upload_type = upload_type

    @staticmethod
    def report_known_error(errors: List[SyncError]):
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

    def upload(self):
        """
        NOTE 此函数运行在其他线程
        上传事件处理
        将收集到的所有上传事件统一触发，上传日志信息
        所有的请求都是网络请求，因此需要异步处理
        对于执行的任务，除了column_type为同步执行，其他为并发执行
        """
        # 根据日志类型进行降重
        upload_tasks_dict = {x: [] for x in self.upload_type}
        # 检查每一个上传结果
        success_tasks_type = []
        # 已知错误列表
        known_errors = []
        # ---------------------------------- 处理upload任务 ----------------------------------

        # 聚合所有的上传任务
        for msg in self.container:
            if msg[0] in self.upload_type:
                upload_tasks_dict[msg[0]].extend(msg[1])

        tasks_key_list = [key for key in upload_tasks_dict if len(upload_tasks_dict[key]) > 0]

        # 同步执行所有的上传任务
        results = [x.value['upload'](upload_tasks_dict[x]) for x in tasks_key_list]
        for index, result in enumerate(results):
            # 如果出现已知问题
            _, e = result
            if isinstance(e, SyncError):
                known_errors.append(e)
                continue
            # 如果出现其他问题，没有办法处理，就直接跳过，但是会有警告
            elif e is not None:
                error = f"{tasks_key_list[index].name} error: {e}, it might be a swanlab bug, data will be lost!"
                swanlog.error(error)
                # continue
                # raise e
            # 标记所有已经成功的任务
            success_tasks_type.append(tasks_key_list[index])

        # ---------------------------------- 最后错误处理 ----------------------------------

        self.container = [(x, upload_tasks_dict[x]) for x in upload_tasks_dict if x not in success_tasks_type]
        self.report_known_error(known_errors)

    def task(self, u: ThreadUtil, *args):
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
            try:
                self.upload()
            except Exception as e:
                swanlog.error(f"upload error: {e}")
            self.lock = False

    def callback(self, u: ThreadUtil, *args):
        """
        回调函数，用于结束时的回调
        NOTE 此函数运行在主线程
        :param u: 线程工具类
        """
        # 如果当前上传任务正在进行，等待上传任务结束
        while self.lock:
            time.sleep(0.1)
        self.container.extend(u.queue.get_all())
        return self.upload()
