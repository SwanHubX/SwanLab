#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:51
@File: http.py
@IDE: pycharm
@Description:
    http会话对象
"""
import httpx
from typing import Optional, Tuple
from datetime import datetime
from .info import LoginInfo, ProjectInfo, ExperimentInfo
from .auth.login import login_by_key
from swanlab.error import NetworkError
from swanlab.package import get_host_api
from swanlab.utils import FONT
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
        self.base_url = get_host_api()
        # 当前cos信息
        self.__cos = None
        # 当前项目信息
        self.__proj = None
        # 当前实验信息
        self.__exp = None

    @property
    def username(self):
        return self.__login_info.username

    @property
    def sid_expired_at(self):
        """
        获取sid的过期时间，字符串格式转时间
        """
        return datetime.strptime(self.__login_info.expired_at, '%Y-%m-%dT%H:%M:%S.%fZ')

    def __create_session(self) -> httpx.AsyncClient:
        """
        创建会话，这将在HTTP类实例化时调用
        """
        session = httpx.AsyncClient(cookies={'sid': self.__login_info.sid})
        # 响应拦截器
        session.event_hooks["response"].append(self.response_interceptor)
        return session

    def __before_request(self):
        # 判断是否已经达到了过期时间
        if (self.sid_expired_at - datetime.utcnow()).total_seconds() < self.REFRESH_TIME:
            # 刷新sid
            self.__login_info = asyncio.run(login_by_key(self.__login_info.api_key))
            self.__session = self.__create_session()

    async def post(self, url: str, data: dict = None) -> dict:
        """
        post请求
        """
        self.__before_request()
        url = self.base_url + url
        resp = await self.__session.post(url, json=data)
        return resp.json()

    def mount_project(self, name):
        async def get_project_info(name: str):
            resp = await http.post(f'/project/{http.username}', data={'name': name})
            return ProjectInfo(resp)

        project: ProjectInfo = asyncio.run(FONT.loading('Getting project...', get_project_info(name)))
        self.__proj = project

    def mount_exp(self, exp_name, colors: Tuple[str, str], description: str = None):
        """
        初始化实验，获取存储信息
        :param exp_name: 所属实验名称
        :param colors: 实验颜色，有两个颜色
        :param description: 实验描述
        """

        async def _():
            """
            创建实验，生成
            :return:
            """
            data = await self.post(f'/project/{self.__login_info.username}/')
            return ExperimentInfo(data)

        self.__exp = asyncio.run(FONT.loading('Creating experiment...', _()))

    @staticmethod
    def response_interceptor(response: httpx.Response):
        """
        捕获所有的http不为2xx的错误
        """
        # 在装饰器中调用被装饰的异步函数
        if response.status_code // 100 != 2:
            raise RuntimeError("http error: {}, reason: {}".format(response.status_code, response.reason_phrase))
        # 返回结果
        return response


http: Optional["HTTP"] = None
"""
一个进程只有一个http请求对象
"""


def create_http(login_info: LoginInfo) -> HTTP:
    """
    创建http请求对象
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
    在使用http对象做请求时，进行统一的错误捕捉
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


def async_refresh_tokens(func):
    """
    捕获错误，刷新http凭证
    """

    async def wrapper(*args, **kwargs):
        return

    return wrapper
