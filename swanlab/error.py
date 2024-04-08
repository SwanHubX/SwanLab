#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:17:14
@File: swanlab/error.py
@IDE: vscode
@Description:
    swanlab全局错误定义
"""
import httpx
from typing import Optional


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


class SyncError(Exception):
    """
    上传错误，作为已知错误捕捉
    """

    def __init__(self, *args):
        super().__init__(*args)
        self.log_level = "error"
        self.message = "sync error"


class ApiError(SyncError):
    """
    api有关的错误，在聚合器中将捕获他们
    """

    def __init__(self, resp: httpx.Response = None, *args):
        super().__init__(*args)
        self.resp: Optional["httpx.Response"] = resp
        self.log_level = "error"
        self.message = 'swanlab api error' if resp is None else 'swanlab api error, status code: {}, reason: {}'.format(
            resp.status_code, resp.reason_phrase
        )


class NetworkError(SyncError):
    """
    请求时网络错误，断网了
    """

    def __init__(self, *args):
        super().__init__(*args)
        self.log_level = "warning"
        self.message = "network error, swanlab will resume uploads when the network improves"
