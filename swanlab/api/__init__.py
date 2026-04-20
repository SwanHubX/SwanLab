"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/20
@description: SwanLab 公共查询 API 入口，面向用户的 OOP 查询接口
"""

from typing import Optional

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.pkg import nrc, scope
from swanlab.sdk.internal.pkg.client import Client
from swanlab.sdk.internal.settings import settings as global_settings

from .experiment import Experiment
from .experiments import Experiments
from .project import Project
from .projects import Projects
from .user import User
from .users import Users
from .workspace import Workspace
from .workspaces import Workspaces


class Api:
    """
    SwanLab 公共查询 API 入口。

    通过独立的 Client 实例与 SwanLab 云端交互，不与 SDK 运行时单例共享。

    用法::

        from swanlab import Api

        api = Api()                              # 自动从 .netrc 读取凭证
        api = Api(api_key="...", host="...")      # 显式传入凭证

        project = api.project("username/project")
        run = api.run("username/project/run_id")
        workspace = api.workspace("username")
        user = api.user()
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        host: Optional[str] = None,
        web_host: Optional[str] = None,
    ):
        """
        初始化 Api 实例。

        认证优先级：显式参数 > Settings（含 .netrc / 环境变量）

        :param api_key: API 密钥，为 None 时从 Settings / .netrc / 环境变量读取
        :param host: API 主机地址，为 None 时从 Settings 读取
        :param web_host: Web 面板地址，为 None 时从 Settings 读取
        """
        api_key, api_host, resolved_web_host = self._resolve_credentials(api_key, host, web_host)
        self._client: Client = Client(api_key=str(api_key), base_url=api_host)
        self._web_host: str = resolved_web_host
        self._api_host = api_host
        self._login_resp = scope.get_context("login_resp")
        self._username: str = self._login_resp["userInfo"]["username"] if self._login_resp else ""

    @staticmethod
    def _resolve_credentials(
        api_key: Optional[str],
        host: Optional[str],
        web_host: Optional[str],
    ) -> tuple:
        """
        按优先级解析凭证：显式参数 > Settings（含 .netrc / 环境变量）。
        返回 (api_key, api_host, web_host)。
        """
        if api_key is None:
            api_key = global_settings.api_key
        if api_key is None:
            raise AuthenticationError("No API key found. Please login with `swanlab login` or pass api_key parameter.")

        api_host: str = nrc.fmt(host) if host is not None else global_settings.api_host
        resolved_web_host: str = nrc.fmt(web_host) if web_host is not None else global_settings.web_host

        return api_key, api_host, resolved_web_host

    # ------------------------------------------------------------------
    #  实体查询方法
    # ------------------------------------------------------------------

    def workspace(self, username: Optional[str] = None) -> Workspace:
        """
        获取工作空间信息，默认为当前登录用户的工作空间。

        :param username: 指定工作空间用户名，为 None 时使用当前登录用户
        """
        if username is None:
            username = self._username
        data = self._client.get(f"/group/{username}").data
        return Workspace(self._client, self._web_host, self._api_host, username=username, data=data)

    def workspaces(self, username: Optional[str] = None) -> Workspaces:
        """
        获取工作空间列表迭代器。

        :param username: 指定用户名，为 None 时使用当前登录用户
        """
        if username is None:
            username = self._username
        return Workspaces(self._client, self._web_host, self._api_host, username=username)

    def project(self, path: str) -> Project:
        """
        获取项目信息。

        :param path: 项目路径，格式为 'username/project-name'
        """
        data = self._client.get(f"/project/{path}").data
        return Project(self._client, self._web_host, self._api_host, path=path, data=data)

    def projects(
        self,
        path: str,
        sort: Optional[str] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ) -> Projects:
        """
        获取工作空间下的项目列表迭代器。

        :param path: 工作空间名称 'username'
        :param sort: 排序方式
        :param search: 搜索关键词
        :param detail: 是否返回详细信息
        """
        return Projects(
            self._client, self._web_host, self._api_host, path=path, sort=sort, search=search, detail=detail
        )

    def run(self, path: str) -> Experiment:
        """
        获取单个实验。

        :param path: 实验路径，格式为 'username/project/run_id'
        """
        parts = path.split("/")
        if len(parts) != 3:
            raise ValueError(f"Invalid path '{path}'. Expected format: 'username/project/run_id'")
        proj_path = path.rsplit("/", 1)[0]
        expid = parts[2]
        data = self._client.get(f"/project/{proj_path}/runs/{expid}").data
        return Experiment(self._client, self._web_host, self._api_host, path=proj_path, data=data)

    def runs(self, path: str, filters: Optional[dict] = None) -> Experiments:
        """
        获取项目下的实验列表迭代器。

        :param path: 项目路径，格式为 'username/project'
        :param filters: 筛选条件
        """
        return Experiments(self._client, self._web_host, self._api_host, path=path, filters=filters)

    def user(self, username: Optional[str] = None) -> User:
        """
        获取用户信息，默认为当前登录用户。

        :param username: 指定用户名
        """
        return User(
            self._client,
            self._web_host,
            self._api_host,
            username=username or self._username,
            login_user=self._username,
        )

    def users(self) -> Users:
        """
        获取用户列表迭代器（私有化部署管理员限定）。
        """
        return Users(self._client, self._web_host, self._api_host, login_user=self._username)


__all__ = ["Api"]
