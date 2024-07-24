#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/19 14:30
@File: utils.py
@IDE: pycharm
@Description:
    任务相关工具函数
"""
from swanlab.package import get_key
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
    def __init__(self, username, src, key, python, name, index):
        """
        :param username: 用户username
        :param key: 用户的api_key
        :param src: 任务zip文件路径
        :param python: 任务的python版本
        :param name: 任务名称
        :param index: 任务入口文件
        """
        self.username = username
        self.src = src
        self.key = key
        self.python = python
        self.name = name
        self.index = index

    def __dict__(self):
        return {
            "username": self.username,
            "src": self.src,
            "index": self.index,
            "python": self.python,
            "conf": {"key": self.key}
        }
