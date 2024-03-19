#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:17:14
@File: swanlab/error.py
@IDE: vscode
@Description:
    swanlab全局错误定义
"""

from .db.error import *


class ValidationError(Exception):
    """验证错误，此时后端验证用户的token或者api key失败

    Parameters
    ----------
    Exception : class
        python内置异常类
    """

    pass


class TokenFileError(Exception):
    """token文件错误，此时token文件不存在或者格式错误（解析失败）

    Parameters
    ----------
    Exception : class
        python内置异常类
    """

    pass


class UnKnownSystemError(Exception):
    """未知系统错误，此时swanlab运行在未知系统上，这个系统不是windows或者类unix系统

    Parameters
    ----------
    Exception : class
        python内置异常类
    """

    pass
