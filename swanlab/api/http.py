#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/7 16:51
@File: http.py
@IDE: pycharm
@Description:
    http会话对象
"""
from typing import Optional, Tuple, Dict, Union, List, AnyStr
from datetime import datetime
from .info import LoginInfo, ProjectInfo, ExperimentInfo
from .auth.login import login_by_key
from .cos import CosClient
from swanlab.data.modules import MediaBuffer
from swanlab.error import NetworkError, ApiError
from swanlab.package import get_host_api, get_package_version
from swankit.log import FONT
from swanlab.log import swanlog
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import json


def decode_response(resp: requests.Response) -> Union[Dict, AnyStr]:
    """
    解码响应，返回信息
    低版本requests库没有JSONDecodeError，所以需要捕获两种异常
    """
    try:
        return resp.json()
    except json.decoder.JSONDecodeError:
        return resp.text
    except requests.JSONDecodeError:
        return resp.text


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
        self.base_url = get_host_api()
        # 当前cos信息
        self.__cos: Optional[CosClient] = None
        # 当前项目信息
        self.__proj: Optional[ProjectInfo] = None
        # 当前实验信息
        self.__exp: Optional[ExperimentInfo] = None
        # 当前进程会话
        self.__session: Optional[requests.Session] = None
        # 当前项目所属的username
        self.__username = login_info.username
        self.__version = get_package_version()
        # 创建会话
        self.__create_session()

    @property
    def groupname(self):
        """
        当前项目所属组名
        """
        return self.__username

    @property
    def username(self):
        """
        当前登录的用户名
        """
        return self.__login_info.username

    @property
    def cos(self):
        return self.__cos

    @property
    def proj_id(self):
        return self.__proj.cuid

    @property
    def projname(self):
        return self.__proj.name

    @property
    def exp_id(self):
        return self.__exp.cuid

    @property
    def expname(self):
        return self.__exp.name

    @property
    def sid_expired_at(self):
        """
        获取sid的过期时间，字符串格式转时间
        """
        return datetime.strptime(self.__login_info.expired_at, "%Y-%m-%dT%H:%M:%S.%fZ")

    def __before_request(self):
        """
        请求前的钩子
        """
        if (self.sid_expired_at - datetime.utcnow()).total_seconds() <= self.REFRESH_TIME:
            # 刷新sid，新建一个会话
            swanlog.debug("Refresh sid...")
            self.__login_info = login_by_key(self.__login_info.api_key, save=False)
            self.__session.headers["cookie"] = f"sid={self.__login_info.sid}"

    def __create_session(self):
        """
        创建会话，这将在HTTP类实例化时调用
        添加了重试策略
        """
        session = requests.Session()
        retry = Retry(
            total=3,
            backoff_factor=0.1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=frozenset(["GET", "POST", "PUT", "DELETE", "PATCH"]),
        )
        adapter = HTTPAdapter(max_retries=retry)
        session.mount("https://", adapter)

        session.headers["swanlab-sdk"] = self.__version
        session.cookies.update({"sid": self.__login_info.sid})

        # 注册响应钩子
        def response_interceptor(response: requests.Response, *args, **kwargs):
            """
            捕获所有的http不为2xx的错误，以ApiError的形式抛出
            """
            if response.status_code // 100 != 2:
                raise ApiError(response, response.status_code, response.reason)

        session.hooks["response"] = response_interceptor

        self.__session = session

    def post(self, url: str, data: dict = None) -> Union[dict, str]:
        """
        post请求
        """
        url = self.base_url + url
        self.__before_request()
        resp = self.__session.post(url, json=data)
        return decode_response(resp)

    def put(self, url: str, data: dict = None) -> Union[dict, str]:
        """
        put请求
        """
        url = self.base_url + url
        self.__before_request()
        resp = self.__session.put(url, json=data)
        return decode_response(resp)

    def get(self, url: str, params: dict = None) -> dict:
        """
        get请求
        """
        url = self.base_url + url
        self.__before_request()
        resp = self.__session.get(url, params=params)
        return decode_response(resp)

    def patch(self, url: str, data: dict = None) -> Union[dict, str]:
        """
        patch请求
        """
        url = self.base_url + url
        self.__before_request()
        resp = self.__session.patch(url, json=data)
        return decode_response(resp)

    def __get_cos(self):
        cos = self.get(f"/project/{self.groupname}/{self.projname}/runs/{self.exp_id}/sts")
        self.__cos = CosClient(cos)

    def upload(self, buffer: MediaBuffer):
        """
        上传文件，需要注意的是file_path应该为unix风格而不是windows风格
        :param buffer: 自定义文件内存对象
        """
        if self.__cos.should_refresh:
            self.__get_cos()
        return self.__cos.upload(buffer)

    def upload_files(self, buffers: List[MediaBuffer]) -> Dict[str, Union[bool, List]]:
        """
        批量上传文件，keys和local_paths的长度应该相等
        :param buffers: 文件内存对象
        :return: 返回上传结果, 包含success_all和detail两个字段，detail为每一个文件的上传结果（通过index索引对应）
        """
        if self.__cos.should_refresh:
            swanlog.debug("Refresh cos...")
            self.__get_cos()
        return self.__cos.upload_files(buffers)

    def mount_project(self, name: str, username: str = None, public: bool = None) -> ProjectInfo:
        """
        创建项目，如果项目已存在，则获取项目信息
        :param name: 项目名称
        :param username: 项目所属的用户名
        :param public: 项目是否公开
        :return: 项目信息
        """
        self.__username = self.__username if username is None else username

        def _():
            try:
                visibility = "PUBLIC" if public else "PRIVATE"
                resp = http.post(f"/project/{self.groupname}", data={"name": name, "visibility": visibility})
            except ApiError as e:
                # 如果为409，表示已经存在，获取项目信息
                if e.resp.status_code == 409:
                    resp = http.get(f"/project/{http.groupname}/{name}")
                elif e.resp.status_code == 404:
                    # 组织/用户不存在
                    raise ValueError(f"Entity `{http.groupname}` not found")
                elif e.resp.status_code == 403:
                    # 权限不足
                    raise ValueError(f"Entity permission denied: " + http.groupname)
                else:
                    raise e
            return ProjectInfo(resp)

        project: ProjectInfo = FONT.loading("Getting project...", _)
        self.__proj = project
        return project

    def mount_exp(self, exp_name, colors: Tuple[str, str], description: str = None):
        """
        初始化实验，获取存储信息
        :param exp_name: 所属实验名称
        :param colors: 实验颜色，有两个颜色
        :param description: 实验描述
        """

        def _():
            """
            先创建实验，后生成cos凭证
            :return:
            """

            data = self.post(
                f"/project/{self.groupname}/{self.__proj.name}/runs",
                {"name": exp_name, "colors": list(colors), "description": description} if description else {
                    "name": exp_name, "colors": list(colors)}
            )
            self.__exp = ExperimentInfo(data)
            # 获取cos信息
            self.__get_cos()

        FONT.loading("Creating experiment...", _)

    def update_state(self, success: bool):
        """
        更新实验状态
        :param success: 实验是否成功
        """

        def _():
            self.put(
                f"/project/{self.groupname}/{self.projname}/runs/{self.exp_id}/state",
                {"state": "FINISHED" if success else "CRASHED", "from": "sdk"},
            )

        FONT.loading("Updating experiment status...", _)


http: Optional["HTTP"] = None
"""
一个进程只有一个http请求对象
"""


def create_http(login_info: LoginInfo) -> HTTP:
    """
    创建http请求对象
    """
    global http
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


def reset_http():
    """
    重置http对象
    """
    global http
    http = None


def sync_error_handler(func):
    """
    在一些接口中我们不希望线程奔溃，而是返回一个错误对象
    """

    def wrapper(*args, **kwargs) -> Tuple[Optional[Union[dict, str]], Optional[Exception]]:
        try:
            # 在装饰器中调用被装饰的异步函数
            result = func(*args, **kwargs)
            return result, None
        except requests.exceptions.Timeout:
            return None, NetworkError()
        except requests.exceptions.ConnectionError:
            return None, NetworkError()
        except Exception as e:
            return None, e

    return wrapper
