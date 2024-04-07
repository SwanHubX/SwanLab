#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:51
@File: http.py
@IDE: pycharm
@Description:
    http会话对象
"""
import requests
from requests.exceptions import RequestException
from typing import Optional
from datetime import datetime
from .auth import LoginInfo
from .auth.login import login_by_key
from swanlab.error import NetworkError
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
"""
一个进程只有一个http请求对象
"""


def init_http(login_info: LoginInfo) -> HTTP:
    """
    初始化http请求对象
    """
    global http
    if http is None:
        http = HTTP(login_info)
    return http


def get_http() -> HTTP:
    """
    创建http请求对象
    :return: http请求对象
    """
    global http
    if http is None:
        raise ValueError("http object is not initialized")
    return http


def async_error_handler(func):
    """
    用于进行统一的错误捕获
    """

    async def wrapper(*args, **kwargs):
        try:
            # 在装饰器中调用被装饰的异步函数
            result = await func(*args, **kwargs)
            return result
        except RequestException:
            return NetworkError()
        except Exception as e:
            return e

    return wrapper
