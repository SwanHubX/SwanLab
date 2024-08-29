#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/19 14:30
@File: utils.py
@IDE: pycharm
@Description:
    任务相关工具函数
"""
from swanlab.package import get_experiment_url
from datetime import datetime, timedelta
import click
import time


def validate_six_char_string(_, __, value):
    if value is None:
        raise click.BadParameter('Parameter is required')
    if not isinstance(value, str):
        raise click.BadParameter('Value must be a string')
    if len(value) != 6:
        raise click.BadParameter('String must be exactly 6 characters long')
    return value


class TaskModel:
    """
    获取到的任务列表模型
    """

    def __init__(self, username: str, task: dict):
        self.cuid = task["cuid"]
        self.username = username
        self.name = task["name"]
        """
        任务名称
        """
        self.python = task["python"]
        """
        任务的python版本
        """
        self.index = task["index"]
        """
        任务的入口文件
        """
        self.project_name = task.get("pName", None)
        """
        项目名称
        """
        self.experiment_id = task.get("eId", None)
        """
        实验ID
        """
        self.created_at = self.fmt_time(task["createdAt"])
        self.started_at = self.fmt_time(task.get("startedAt", None))
        self.finished_at = self.fmt_time(task.get("finishedAt", None))
        self.status = task["status"]
        self.msg = task.get("msg", None)
        self.combo = task["combo"]

    @property
    def url(self):
        if self.project_name is None or self.experiment_id is None:
            return None
        return get_experiment_url(self.username, self.project_name, self.experiment_id)

    @staticmethod
    def fmt_time(date: str = None):
        if date is None:
            return None
        date = date.replace("Z", "+00:00")
        # 获取当前计算机时区的时差（以秒为单位）
        local_time_offset = time.altzone if time.localtime().tm_isdst else time.timezone
        local_time_offset = timedelta(seconds=-local_time_offset)
        local_time = datetime.fromisoformat(date) + local_time_offset
        # 将 UTC 时间转换为本地时间
        return local_time.strftime("%Y-%m-%d %H:%M:%S")
