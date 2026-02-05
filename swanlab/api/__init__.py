"""
@author: Zhou QiYang
@file: __init__.py
@time: 2026/1/5 17:58
@description: SwanLab OpenAPI包
"""

from typing import Optional, List, Dict

from swanlab.core_python import auth, Client
from swanlab.core_python.api.experiment import get_single_experiment, get_project_experiments
from swanlab.error import KeyFileError
from swanlab.log import swanlog
from swanlab.package import HostFormatter, get_key
from .deprecated import OpenApi
from .experiment import Experiment
from .experiments import Experiments
from .project import Project
from .projects import Projects
from .user import User
from .users import Users
from .utils import self_hosted
from .workspace import Workspace
from .workspaces import Workspaces
from ..core_python.api.project import get_project_info
from ..core_python.api.user import get_workspace_info


class Api:
    def __init__(self, api_key: Optional[str] = None, host: Optional[str] = None, web_host: Optional[str] = None):
        """
        初始化 OpenApi 实例，用户需提前登录，或者提供API密钥
        :param api_key: API 密钥，可选
        :param host: API 主机地址，可选
        :param web_host: Web 主机地址，可选
        """
        if host or web_host:
            HostFormatter(host, web_host)()
        if api_key:
            swanlog.debug("Using API key", api_key)
        else:
            swanlog.debug("Using existing key")
            try:
                api_key = get_key()
            except KeyFileError as e:
                swanlog.error("To use SwanLab OpenAPI, please login first.")
                raise RuntimeError("Not logged in.") from e

        self._login_info = auth.code_login(api_key, save_key=False)
        # 一个OpenApi对应一个client，可创建多个api获取从不同的client获取不同账号下的实验信息
        self._client: Client = Client(self._login_info)
        self._web_host = self._login_info.web_host
        self._login_user = self._login_info.username

    @self_hosted("root")
    def users(self) -> Users:
        """
        超级管理员获取所有用户
        :return: User 实例，可对当前/指定用户进行操作
        """
        return Users(self._client, login_user=self._login_user)

    def user(self, username: str = None) -> User:
        """
        获取用户实例，用于操作用户相关信息
        :param username: 指定用户名，如果为 None，则返回当前登录用户
        :return: User 实例，可对当前/指定用户进行操作
        """
        return User(client=self._client, login_user=self._login_user, username=username)

    def workspaces(
        self,
        username: str = None,
    ):
        """
        获取当前登录用户的工作空间迭代器
        当username为其他用户时，可以作为visitor访问其工作空间
        """
        if username is None:
            username = self._login_user
        return Workspaces(self._client, username=username)

    def workspace(
        self,
        username: str = None,
    ):
        """
        获取当前登录用户的工作空间
        """
        if username is None:
            username = self._login_user
        data = get_workspace_info(self._client, path=username)
        return Workspace(self._client, data=data)

    def projects(
        self,
        path: str,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ) -> Projects:
        """
        获取指定工作空间（组织）下的所有项目信息
        :param path: 工作空间（组织）名称 'username'
        :param sort: 排序方式，可选
        :param search: 搜索关键词，可选
        :param detail: 是否返回详细信息，可选
        :return: Projects 实例，可遍历获取项目信息
        """
        return Projects(
            self._client,
            web_host=self._web_host,
            path=path,
            sort=sort,
            search=search,
            detail=detail,
        )

    def project(
        self,
        path: str,
    ) -> Project:
        """
        获取指定工作空间（组织）下的指定项目信息
        :param path: 项目路径 'username/project'
        :return: Project 实例，单个项目的信息
        """
        data = get_project_info(self._client, path=path)
        return Project(self._client, web_host=self._web_host, data=data, login_info=self._login_info)

    def runs(self, path: str, filters: Dict[str, object] = None) -> Experiments:
        """
        获取指定项目下的所有实验信息
        :param path: 项目路径，格式为 'username/project'
        :return: Experiments 实例，可遍历获取实验信息
        :param filters: 筛选实验的条件，可选
        """
        return Experiments(self._client, path=path, login_info=self._login_info, filters=filters)

    def run(
        self,
        path: str,
    ) -> Experiment:
        """
        获取指定实验的信息
        :param path: 实验路径，格式为 'username/project/run_id'
        :return: Experiment 实例，包含实验信息
        """
        if len(path.split('/')) != 3:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project/run_id'")
        data = get_single_experiment(self._client, path=path)
        proj_path = path.rsplit('/', 1)[0]
        return Experiment(
            self._client,
            data=data,
            path=proj_path,
            web_host=self._web_host,
            login_user=self._login_user,
            line_count=1,
        )


__all__ = ["Api", "OpenApi"]