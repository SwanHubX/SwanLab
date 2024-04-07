#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 21:52
@File: __init__.py.py
@IDE: pycharm
@Description:
    API模块，封装api请求接口
"""
from .auth import LoginInfo
import requests
from typing import Optional
from .auth.login import login_request, terminal_login, code_login, login_by_key
from swanlab.error import ValidationError
from datetime import datetime
import asyncio


class HTTP:
    """
    封装请求函数，添加get、post、put、delete方法
    """
    REFRESH_TIME = 60 * 60 * 24 * 7  # 7天
    """
    刷新时间，单位秒，如果sid过期时间减去当前时间小于这个时间，就刷新sid
    """

    def __init__(self, login_info: LoginInfo):
        """
        初始化会话
        """
        self.__login_info = login_info
        self.__session = self.__create_session()

    def expired_at(self):
        """
        获取sid的过期时间，字符串格式转时间
        """
        return datetime.strptime(self.__login_info.expired_at, '%Y-%m-%dT%H:%M:%S.%fZ')

    def __create_session(self) -> requests.Session:
        """
        创建会话，这将在HTTP类实例化时调用
        """
        req = requests.Session()
        req.cookies.update({'sid': self.__login_info.sid})
        return req

    def __before_request(self):
        # 判断是否已经达到了过期时间
        if (self.expired_at() - datetime.utcnow()).total_seconds() < self.REFRESH_TIME:
            # 刷新sid
            self.__login_info = asyncio.run(login_by_key(self.__login_info.api_key))
            self.__session = self.__create_session()

    async def get(self, url: str, **kwargs) -> requests.Response:
        """
        get请求
        """
        self.__before_request()
        return self.__session.get(url, **kwargs)

    async def post(self, url: str, **kwargs) -> requests.Response:
        """
        post请求
        """
        self.__before_request()
        return self.__session.post(url, **kwargs)

    async def put(self, url: str, **kwargs) -> requests.Response:
        """
        put请求
        """
        self.__before_request()
        return self.__session.put(url, **kwargs)

    async def delete(self, url: str, **kwargs) -> requests.Response:
        """
        delete请求
        """
        self.__before_request()
        return self.__session.delete(url, **kwargs)


http: Optional["HTTP"] = None


def create_http(login_info: LoginInfo) -> HTTP:
    """
    创建http请求对象
    :return: http请求对象
    """
    global http
    http = HTTP(login_info)
    return http


__all__ = ["LoginInfo", "code_login", 'terminal_login', "HTTP", "create_http"]
