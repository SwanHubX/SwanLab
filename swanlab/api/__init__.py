#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 21:52
@File: __init__.py.py
@IDE: pycharm
@Description:
    API模块，封装api请求接口
"""
from .auth.login import terminal_login, code_login
from .http import create_http, get_http
from .info import *


def is_login() -> bool:
    """判断是否登录(拥有http对象)"""
    try:
        return get_http() is not None
    except ValueError:
        return False


__all__ = [
    "LoginInfo",
    'ExperimentInfo',
    'ProjectInfo',
    "code_login",
    'terminal_login',
    'create_http',
    'get_http',
    'is_login',
]
