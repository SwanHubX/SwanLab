"""
@author: cunyue
@file: __init__.py
@time: 2025/6/16 13:29
@description: swanlab 客户端，负责发送 http 请求
"""

import json
from datetime import datetime, timezone
from typing import Optional, Tuple, Dict, Union, List, AnyStr

import requests
from urllib3.exceptions import (
    MaxRetryError,
    TimeoutError,
    NewConnectionError,
    ConnectionError,
    ReadTimeoutError,
    ConnectTimeoutError,
)

from swanlab.error import NetworkError, ApiError
from swanlab.log import swanlog
from swanlab.package import get_package_version
from .model import ProjectInfo, ExperimentInfo
from .session import create_session
from .. import auth
from ...env import utc_time


def decode_response(resp: requests.Response) -> Union[Dict, AnyStr, List]:
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


class Client:
    """
    封装请求函数，添加get、post、put、delete方法
    会自动刷新会话信息
    """

    REFRESH_TIME = 60 * 60 * 24 * 7  # 7天
    """
    刷新时间，单位秒，如果sid过期时间减去当前时间小于这个时间，就刷新sid
    """

    def __init__(self, login_info: auth.LoginInfo):
        self.__login_info = login_info
        # 当前会话
        self.__session: Optional[requests.Session] = None
        self.__version = get_package_version()
        self.__create_session()

        # 标识当前实验会话（flagId）是否被其他进程顶掉
        self.pending = False
        # 当前项目所属的username
        self.__groupname = login_info.username
        # 当前项目信息
        self.__proj: Optional[ProjectInfo] = None
        # 当前实验信息
        self.__exp: Optional[ExperimentInfo] = None

    # ---------------------------------- 一些辅助属性 ----------------------------------
    @property
    def exp(self) -> ExperimentInfo:
        assert self.__exp is not None, "Experiment not mounted, please call mount_exp() first"
        return self.__exp

    @property
    def proj(self) -> ProjectInfo:
        assert self.__proj is not None, "Project not mounted, please call mount_project() first"
        return self.__proj

    @property
    def api_key(self):
        return self.__login_info.api_key

    @property
    def groupname(self):
        """
        当前项目所属组名
        """
        return self.__groupname

    @property
    def username(self):
        """
        当前登录的用户名
        """
        return self.__login_info.username

    @property
    def projname(self):
        return self.__proj.name

    @property
    def history_exp_count(self):
        return self.__proj.history_exp_count

    @property
    def exp_id(self):
        return self.__exp.cuid

    @property
    def expname(self):
        return self.__exp.name

    @property
    def web_proj_url(self):
        return f"{self.__login_info.web_host}/@{self.groupname}/{self.projname}"

    @property
    def web_exp_url(self):
        return f"{self.web_proj_url}/runs/{self.exp_id}"

    # ---------------------------------- http方法 ----------------------------------

    def __before_request(self):
        """
        请求前的钩子
        """
        sid_expired_at = datetime.strptime(self.__login_info.expired_at, "%Y-%m-%dT%H:%M:%S.%fZ").replace(
            tzinfo=timezone.utc
        )
        if (sid_expired_at - utc_time()).total_seconds() <= self.REFRESH_TIME:
            # 刷新sid，新建一个会话
            swanlog.debug("Refresh sid...")
            self.__login_info = auth.login_by_key(self.__login_info.api_key, save=False)
            self.__session.headers["cookie"] = f"sid={self.__login_info.sid}"

        # 携带实验会话Id
        if self.__exp is not None:
            self.__session.headers["flagId"] = self.__exp.flag_id

    def __create_session(self):
        """
        创建会话，这将在HTTP类实例化时调用
        添加了重试策略
        """
        session = create_session()
        session.headers["swanlab-sdk"] = self.__version
        session.cookies.update({"sid": self.__login_info.sid})

        # 注册响应钩子
        def response_interceptor(response: requests.Response, *args, **kwargs):
            """
            捕获所有的http不为2xx的错误，以ApiError的形式抛出
            """
            # 1. 日志打印
            swanlog.debug(
                f"HTTP Request: {response.request.method.upper()} {response.url} | "
                f"Response Status: {response.status_code} | "
                f"Body: {decode_response(response)}"
            )
            # 2. 如果状态码不为2xx，抛出异常
            if response.status_code // 100 != 2:
                traceid = f"Trace id: {response.headers.get('traceid')}"
                request = f"{response.request.method.upper()} {response.url}"
                resp = f"{response.status_code} {response.reason}"
                raise ApiError(response, traceid, request, resp)

        session.hooks["response"] = response_interceptor

        self.__session = session

    def post(self, url: str, data: Union[dict, list] = None):
        """
        post请求
        """
        url = self.__login_info.api_host + url
        self.__before_request()
        resp = self.__session.post(url, json=data)
        return decode_response(resp), resp

    def put(self, url: str, data: dict = None):
        """
        put请求
        """
        url = self.__login_info.api_host + url
        self.__before_request()
        resp = self.__session.put(url, json=data)
        return decode_response(resp), resp

    def get(self, url: str, params: dict = None):
        """
        get请求
        """
        url = self.__login_info.api_host + url
        self.__before_request()
        resp = self.__session.get(url, params=params)
        return decode_response(resp), resp

    def patch(self, url: str, data: dict = None):
        """
        patch请求
        """
        url = self.__login_info.api_host + url
        self.__before_request()
        resp = self.__session.patch(url, json=data)
        return decode_response(resp), resp

    # ---------------------------------- 训练相关接口 ----------------------------------

    def mount_project(self, name: str, username: str = None, public: bool = None):
        """
        创建项目，如果项目已存在，则获取项目信息
        :param name: 项目名称
        :param username: 项目所属的用户名
        :param public: 项目是否公开
        :return: 项目信息
        """
        try:
            data = {"name": name}
            if username is not None:
                data["username"] = username
            if public is not None:
                data["visibility"] = "PUBLIC" if public else "PRIVATE"
            resp_data, _ = self.post(f"/project", data=data)
        except ApiError as e:
            if e.resp.status_code == 409:
                # 项目已经存在，从对象中解析信息
                resp_data = decode_response(e.resp)
            elif e.resp.status_code == 404 and e.resp.reason == "Not Found":
                # 早期（主要是私有化）swanlab 后端没有 /project 接口，需要使用 /project/{username} 接口，此时没有默认空间的特性
                self.__groupname = self.__groupname if username is None else username
                try:
                    visibility = "PUBLIC" if public else "PRIVATE"
                    resp_data, _ = self.post(
                        f"/project/{self.groupname}", data={"name": name, "visibility": visibility}
                    )
                except ApiError as e:
                    # 如果为409，表示已经存在，获取项目信息
                    if e.resp.status_code == 409:
                        resp_data, _ = self.get(f"/project/{self.groupname}/{name}")
                    elif e.resp.status_code == 404:
                        # 组织/用户不存在
                        raise ValueError(f"Space `{self.groupname}` not found")
                    elif e.resp.status_code == 403:
                        # 权限不足
                        raise ValueError(f"Space permission denied: " + self.groupname)
                    else:
                        raise e
            else:
                # 此接口为后端处理，sdk 在理论上不会出现其他错误，因此不需要处理其他错误
                raise e
        # 设置当前项目所属的用户名
        self.__groupname = resp_data['username']
        # 获取详细信息
        resp_data_info, _ = self.get(f"/project/{self.groupname}/{name}")
        self.__proj = ProjectInfo(resp_data_info)

    def mount_exp(
        self,
        exp_name,
        colors: Tuple[str, str],
        description: str = None,
        job_type: str = None,
        group: str = None,
        tags: List[str] = None,
        created_at: str = None,
        cuid: str = None,
        must_exist: bool = False,
    ) -> bool:
        """
        初始化实验，获取存储信息
        :param exp_name: 所属实验名称
        :param colors: 实验颜色，有两个颜色
        :param description: 实验描述
        :param job_type: 任务类型
        :param group: 实验组
        :param tags: 实验标签
        :param created_at: 实验创建时间，格式为 ISO 8601
        :param cuid: 实验的唯一标识符，如果不提供则由后端生成
        :param must_exist: 如果 cuid 被传递，是否限制实验必须存在

        :raises RuntimeError: 如果实验不存在且must_exist为True
        :raises NotImplementedError: 如果项目未挂载

        :return: 返回实验为新建的还是更新的，为 True 时为新建实验
        """
        assert self.__proj is not None, "Project not mounted, please call mount_project() first"
        if must_exist:
            assert cuid is not None, "cuid must be provided when must_exist is True"
            try:
                self.get(f"/project/{self.groupname}/{self.__proj.name}/runs/{cuid}")
            except ApiError as e:
                if e.resp.status_code == 404 and e.resp.reason == "Not Found":
                    raise RuntimeError(f"Experiment {cuid} does not exist in project {self.projname}")

        labels = [{"name": tag} for tag in tags] if tags else []
        post_data = {
            "name": exp_name,
            "description": description,
            "createdAt": created_at,
            "colors": list(colors),
            "labels": labels if len(labels) else None,
            "job": job_type,
            "cluster": group,
            "cuid": cuid,
        }
        post_data = {k: v for k, v in post_data.items() if v is not None}  # 移除值为None的键

        # 这部分错误将不会被上层捕获，直接抛出异常
        try:
            data, resp = self.post(f"/project/{self.groupname}/{self.__proj.name}/experiment", post_data)
        except ApiError as e:
            if e.resp.status_code == 400 and e.resp.reason == "Bad Request":
                # 指定的 cuid 对应的实验是克隆实验
                raise ValueError(
                    f"Experiment with CUID {cuid} is a cloned experiment (cloned experiments cannot be resumed).",
                )
            elif e.resp.status_code == 403 and e.resp.reason == "Forbidden":
                # 权限不足
                raise ValueError(f"Project permission denied: {self.projname}")
            elif e.resp.status_code == 404 and e.resp.reason == "Not Found":
                # 传入的项目不存在
                raise ValueError(f"Project {self.projname} not found")
            elif e.resp.status_code == 404 and e.resp.reason == "Disabled Resource":
                # 传入的实验被删除
                raise ValueError(f"Experiment {cuid} has been deleted")
            elif e.resp.status_code == 409 and e.resp.reason == "Conflict":
                # 传入 cuid 但是实验不属于当前项目
                raise ValueError(f"Experiment with CUID {cuid} does not belong to project {self.projname}")
            raise e
        # 200代表实验已存在，开启更新模式
        # 201代表实验不存在，新建实验
        new = resp.status_code == 201
        # 这部分信息暂时没有用到
        self.__exp = ExperimentInfo(data)
        # 重置挂起状态
        self.pending = False
        return new


client: Optional["Client"] = None
"""
一个进程只有一个客户端对象
"""


def create_client(login_info: auth.LoginInfo) -> Client:
    """
    创建客户端对象
    """
    global client
    client = Client(login_info)
    return client


def get_client() -> Client:
    """
    获取客户端对象
    :return: client
    """
    global client
    if client is None:
        raise ValueError("client object is not initialized")
    return client


def reset_client():
    """
    重置client对象
    """
    global client
    client = None


def safe_request(func):
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
        # Catch urllib3 specific errors
        except (
            MaxRetryError,
            TimeoutError,
            NewConnectionError,
            ConnectionError,
            ReadTimeoutError,
            ConnectTimeoutError,
        ):
            return None, NetworkError()
        except Exception as e:
            return None, e

    return wrapper


__all__ = [
    "get_client",
    "reset_client",
    "create_session",
    "create_client",
    "safe_request",
    "decode_response",
    "Client",
]
