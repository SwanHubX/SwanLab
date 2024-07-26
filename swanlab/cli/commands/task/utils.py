#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/19 14:30
@File: utils.py
@IDE: pycharm
@Description:
    任务相关工具函数
"""
from swanlab.package import get_key, get_experiment_url
from swanlab.api import terminal_login, create_http, LoginInfo
from swanlab.error import KeyFileError


def login_init_sid() -> LoginInfo:
    key = None
    try:
        key = get_key()
    except KeyFileError:
        pass
    login_info = terminal_login(key)
    create_http(login_info)
    return login_info


class TaskModel:
    """
    获取到的任务列表模型
    """

    def __init__(self, username: str, task: dict, ):
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
        self.project_name = task.get("pName", None)
        """
        项目名称
        """
        self.experiment_id = task.get("eId", None)
        """
        实验ID
        """
        self.created_at = task["createdAt"]
        self.started_at = task.get("startedAt", None)
        self.finished_at = task.get("finishedAt", None)
        self.status = task["status"]
        self.msg = task.get("msg", None)

    @property
    def url(self):
        if self.project_name is None or self.experiment_id is None:
            return None
        return get_experiment_url(self.username, self.project_name, self.experiment_id)
