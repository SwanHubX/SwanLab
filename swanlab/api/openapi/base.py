#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:36
@File: base.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI API基类
"""
from functools import wraps

from swanlab.api import LoginInfo
from swanlab.api.http import HTTP
from swanlab.error import ApiError


def handle_api_error(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except ApiError as e:
            return {"code": e.resp.status_code, "message": e.message}

    return wrapper


class ApiHTTP:
    def __init__(self, login_info: LoginInfo):
        self.__http: HTTP = HTTP(login_info)

    @property
    def username(self):
        """
        当前登录的用户名
        """
        return self.http.username

    @property
    def http(self) -> HTTP:
        """
        当前使用的HTTP对象
        """
        return self.__http

    @handle_api_error
    def get(self, url: str, params: dict = None):
        return self.http.get(url=url, params=params)

    @handle_api_error
    def post(self, url: str, data: dict = None):
        return self.http.post(url=url, data=data)


class ApiBase:
    def __init__(self, http: ApiHTTP):
        self.http: ApiHTTP = http
