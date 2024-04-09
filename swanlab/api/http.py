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
from typing import Optional, Tuple, Dict, Union, List
from datetime import datetime
from .info import LoginInfo, ProjectInfo, ExperimentInfo
from .auth.login import login_by_key
from .cos import CosClient
from swanlab.error import NetworkError, ApiError
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
        self.__cos: Optional["CosClient"] = None
        # 当前项目信息
        self.__proj: Optional["ProjectInfo"] = None
        # 当前实验信息
        self.__exp: Optional["ExperimentInfo"] = None

    @property
    def username(self):
        return self.__login_info.username

    @property
    def proj_id(self):
        return self.__proj.cuid

    @property
    def exp_id(self):
        return self.__exp.cuid

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

        # 注册请求前的钩子
        async def request_interceptor(request: httpx.Request):
            # 判断是否已经达到了过期时间
            if (self.sid_expired_at - datetime.utcnow()).total_seconds() <= self.REFRESH_TIME:
                # 刷新sid，新建一个会话
                self.__login_info = await login_by_key(self.__login_info.api_key)
                self.__session = self.__create_session()
                # 更新当前请求的cookie
                request.headers['cookie'] = f'sid={self.__login_info.sid}'

        session.event_hooks['request'].append(request_interceptor)

        # 注册响应钩子
        async def response_interceptor(response: httpx.Response):
            """
            捕获所有的http不为2xx的错误，以ApiError的形式抛出
            """
            # 如果是
            if response.status_code // 100 != 2:
                raise ApiError(response)

        session.event_hooks['response'].append(response_interceptor)

        return session

    async def post(self, url: str, data: dict = None) -> dict:
        """
        post请求
        """
        url = self.base_url + url
        resp = await self.__session.post(url, json=data)
        return resp.json()

    async def get(self, url: str, params: dict = None) -> dict:
        """
        get请求
        """
        url = self.base_url + url
        resp = await self.__session.get(url, params=params)
        return resp.json()

    async def __get_cos(self):
        cos = await self.get(f'/project/{self.__login_info.username}/{self.__proj.name}/runs/{self.__exp.cuid}/sts')
        self.__cos = CosClient(cos)

    async def upload(self, key: str, local_path):
        """
        上传文件，需要注意的是file_path应该为unix风格而不是windows风格
        开头不能有/，即使有也会被去掉
        :param key: 上传到cos的文件名称
        :param local_path: 本地文件路径，一般用绝对路径
        """
        if key.startswith('/'):
            key = key[1:]
        if self.__cos.should_refresh:
            await self.__get_cos()
        return self.__cos.upload(key, local_path)

    async def upload_files(self, keys: list, local_paths: list) -> Dict[str, Union[bool, List]]:
        """
        批量上传文件，keys和local_paths的长度应该相等
        :param keys: 上传到cos
        :param local_paths: 本地文件路径，需用绝对路径
        :return: 返回上传结果, 包含success_all和detail两个字段，detail为每一个文件的上传结果（通过index索引对应）
        """
        if self.__cos.should_refresh:
            await self.__get_cos()
        keys = [key[1:] if key.startswith('/') else key for key in keys]
        return self.__cos.upload_files(keys, local_paths)

    def mount_project(self, name: str):
        async def _():
            try:
                resp = await http.post(f'/project/{http.username}', data={'name': name})
            except ApiError as e:
                # 如果为409，表示已经存在，获取项目信息
                if e.resp.status_code == 409:
                    resp = await http.get(f'/project/{http.username}/{name}')
                else:
                    raise e
            return ProjectInfo(resp)

        project: ProjectInfo = asyncio.run(FONT.loading('Getting project...', _()))
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
            data = await self.post(
                f'/project/{self.__login_info.username}/{self.__proj.name}/runs',
                {"name": exp_name, "color": list(colors), "description": description}
            )
            self.__exp = ExperimentInfo(data)
            # 获取cos信息
            await self.__get_cos()

        asyncio.run(FONT.loading('Creating experiment...', _()))


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
    在一些接口中我们不希望线程奔溃，而是返回一个错误对象
    """

    async def wrapper(*args, **kwargs):
        try:
            # 在装饰器中调用被装饰的异步函数
            result = await func(*args, **kwargs)
            return result
        except httpx.NetworkError:
            return NetworkError()
        except Exception as e:
            return e

    return wrapper
