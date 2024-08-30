#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:17:14
@File: swanlab/error.py
@IDE: vscode
@Description:
    swanlab全局错误定义
"""
import requests


class APIKeyFormatError(Exception):
    """错误的api key格式，此类承担解析api key的任务
    错误：类型错误、长度错误、字符错误
    """

    API_KEY_LENGTH = 21

    @classmethod
    def check(cls, api_key: str):
        """
        判断是否抛出异常，先检查长度，再检查类型
        """
        if not isinstance(api_key, str):
            raise cls("Api key must be a string")
        for c in api_key:
            # 0-9, A-Z, a-z
            if ord(c) not in range(48, 58) and ord(c) not in range(65, 91) and ord(c) not in range(97, 123):
                raise cls("Invalid character in api key: {}".format(repr(c)))
        if len(api_key) != 21:
            raise cls("Api key length must be 21 characters long, yours was {}".format(len(api_key)))
        return api_key


class ValidationError(Exception):
    """验证错误，此时后端验证用户的token或者api key失败"""

    pass


class KeyFileError(Exception):
    """key存储的文件错误，此时key文件不存在或者格式错误（解析失败）"""

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

    def __init__(self, resp: requests.Response = None, *args):
        super().__init__(*args)
        self.resp = resp
        self.log_level = "error"
        self.message = (
            'swanlab api error'
            if resp is None
            else 'swanlab api error, status code: {}, reason: {}'.format(resp.status_code, resp.reason)
        )


class NetworkError(SyncError):
    """
    请求时网络错误，断网了
    """

    def __init__(self, *args):
        super().__init__(*args)
        self.log_level = "warning"
        self.message = "network error, swanlab will resume uploads when the network improves"


class DataTypeError(Exception):
    """数据类型错误，此时数据类型不符合预期"""

    def __init__(self, expected: str, got: str):
        """
        :param expected: 期望的数据类型
        :param got: 实际的数据类型
        """
        super().__init__(f"expected: {expected}, got: {got}")
        self.expected = expected
        self.got = got
