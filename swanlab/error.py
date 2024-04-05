#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:17:14
@File: swanlab/error.py
@IDE: vscode
@Description:
    swanlab全局错误定义
"""


class ValidationError(Exception):
    """验证错误，此时后端验证用户的token或者api key失败
    """

    pass


class KeyFileError(Exception):
    """key存储的文件错误，此时key文件不存在或者格式错误（解析失败）
    """

    pass


class NotLoginError(Exception):
    """未登录错误，此时用户未登录
    """
    pass


class UnKnownSystemError(Exception):
    """未知系统错误，此时swanlab运行在未知系统上，这个系统不是windows或者类unix系统
    """
    pass


class UpLoadError(Exception):
    """
    日志上传有关的错误，在聚合器中将捕获他们
    """
    pass


class NetworkError(UpLoadError):
    """
    请求时网络错误，断网了
    """
    pass
