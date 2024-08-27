#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/25 14:33
@File: utils.py
@IDE: pycharm
@Description:
    一些工具函数
"""
from typing import Optional
from swanlab.log import swanlog
from swanlab.package import get_key
from swanlab.api import terminal_login, create_http, LoginInfo, get_http
from swanlab.error import KeyFileError, ApiError
import sys


class UseTaskHttp:
    """
    主要用于检测http响应是否为3xx字段，如果是则要求用户更新版本
    使用此类之前需要先调用login_init_sid()函数完成全局http对象的初始化
    """

    def __init__(self):
        self.http = get_http()

    def __enter__(self):
        return self.http

    def __exit__(self, exc_type, exc_val: Optional[ApiError], exc_tb):
        if exc_type is ApiError:
            # api已过期，需要更新swanlab版本
            if exc_val.resp.status_code // 100 == 3:
                swanlog.info("SwanLab in your environment is outdated. Upgrade: `pip install -U swanlab`")
                sys.exit(3)
        return False


def login_init_sid() -> LoginInfo:
    key = None
    try:
        key = get_key()
    except KeyFileError:
        pass
    login_info = terminal_login(key)
    create_http(login_info)
    return login_info
