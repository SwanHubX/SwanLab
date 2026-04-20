"""
@author: Nexisato
@file: __init__.py
@time: 2026/4/20
@description: SwanLab 公共查询 API 入口，面向用户的 OOP 查询接口
"""

from typing import Optional

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.pkg import nrc, scope
from swanlab.sdk.internal.pkg.client import Client
from swanlab.sdk.internal.settings import settings as global_settings
from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

from .experiment import Experiment
from .project import Project
from .user import User
from .workspace import Workspace


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

        认证优先级与 swanlab.login 对齐：显式参数 > Settings（含 .netrc / 环境变量）

        :param api_key: API 密钥，为 None 时从 Settings / .netrc / 环境变量读取
        :param host: API 主机地址，为 None 时从 Settings 读取
        :param web_host: Web 面板地址，为 None 时从 Settings 读取
        """
        api_key, api_host, web_host = self._resolve_credentials(api_key, host, web_host)
        # 创建独立 Client 实例，不走 core_python/client 模块级单例
        self._client: Client = Client(api_key=api_key, base_url=api_host)
        self._web_host: str = web_host
        self._login_resp: Optional[LoginResponse] = scope.get_context("login_resp")
        self._username: str = self._login_resp["userInfo"]["username"] if self._login_resp else ""

    @staticmethod
    def _resolve_credentials(
        api_key: Optional[str],
        host: Optional[str],
        web_host: Optional[str],
    ) -> tuple[str, str, str]:
        """
        按优先级解析凭证，与 swanlab.login 对齐：显式参数 > Settings（含 .netrc / 环境变量）。
        返回 (api_key, api_host, web_host)。
        """
        # api_key：显式传入 > global_settings（已含 .netrc > 环境变量 > 默认值的优先链）
        if api_key is None:
            api_key = global_settings.api_key
        if api_key is None:
            raise AuthenticationError("No API key found. Please login with `swanlab login` or pass api_key parameter.")

        # host：显式传入 > global_settings
        api_host: str = nrc.fmt(host) if host is not None else global_settings.api_host

        # web_host：显式传入 > global_settings
        resolved_web_host: str
        if web_host is not None:
            resolved_web_host = nrc.fmt(web_host)
        else:
            resolved_web_host = global_settings.web_host

        return api_key, api_host, resolved_web_host

    # ------------------------------------------------------------------
    #  实体查询方法
    # ------------------------------------------------------------------

    def workspace(self, username: Optional[str] = None) -> Workspace:
        """获取工作空间信息，默认为当前登录用户的工作空间。"""
        if username is None:
            username = self._username
        return Workspace(self._client, username=username, web_host=self._web_host)

    def project(self, path: str) -> Project:
        """
        获取项目信息。

        :param path: 项目路径，格式为 'username/project-name'
        """
        return Project(self._client, path=path, web_host=self._web_host)

    def user(self, username: Optional[str] = None) -> User:
        """
        获取用户信息，默认为当前登录用户。

        :param username: 指定用户名
        """
        if username is None:
            username = self._username
        return User(self._client, username=username)

    def run(self, path: str) -> Experiment:
        """
        获取单个实验。

        :param path: 实验路径，格式为 'username/project/run_id'
        """
        parts = path.split("/")
        if len(parts) != 3:
            raise ValueError(f"Invalid path '{path}'. Expected format: 'username/project/run_id'")
        return Experiment(self._client, path=path, web_host=self._web_host)


__all__ = ["Api"]
