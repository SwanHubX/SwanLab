#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/4/29 9:40
@File: __init__.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI模块
"""

from swanlab.api import create_http, is_login, code_login
from swanlab.error import KeyFileError, ValidationError
from swanlab.log import swanlog
from swanlab.package import get_key

def pre_login():
    """检查登录状态, 若未登录，则尝试登录"""
    if not is_login():
        try:
            key = get_key()
        except KeyFileError as e:
            swanlog.error("To use SwanLab OpenAPI, please login first.")
            raise RuntimeError("Not logged in.") from e
        login_info = code_login(key, False)
        if login_info.is_fail:
            raise ValidationError("Login failed: " + str(login_info))
        create_http(login_info)

pre_login()

from .group import (
    list_workspaces
)
