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

from .base import ApiClientContext, BaseEntity
from .experiment import Experiment, Experiments
from .project import Project, Projects
from .typings.common import ApiResponseType, PaginatedQuery
from .user import User
from .workspace import Workspace, Workspaces


class Api(BaseEntity):
    """
    SwanLab 公共查询 API 入口。

    通过独立的 Client 实例与 SwanLab 云端交互，与 SDK 运行时客户端完全隔离。
    继承 BaseEntity 以复用 _get/_post/_put/_delete/_paginate 等安全 HTTP 方法。

    用法::

        from swanlab import Api

        api = Api()                              # 自动从 .netrc 读取凭证
        api = Api(api_key="...", host="...")      # 显式传入凭证

        resp = api.project("username/project")
        if resp.ok:
            project = resp.data
            print(project.name)
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        host: Optional[str] = None,
        web_host: Optional[str] = None,
    ) -> None:
        """
        初始化 Api 实例。

        认证优先级：
        1. 显式参数 (api_key / host / web_host)
        2. scope 登录态（进程内已调用 swanlab.login 时可用）
        3. Settings（含 .netrc / 环境变量）

        始终创建独立的 Client 实例，与 SDK 运行时单例互不干扰。

        :param api_key: API 密钥，为 None 时从 Settings / .netrc / 环境变量读取
        :param host: API 主机地址，为 None 时从 Settings 读取
        :param web_host: Web 面板地址，为 None 时从 Settings 读取
        """
        # 优先从 scope 获取已有登录态（如进程内已调用 swanlab.login），直接复用凭证
        login_resp = scope.get_context("login_resp")
        api_key, api_host, web_host = self._resolve_credentials(api_key, host, web_host)
        _client = Client(api_key=str(api_key), base_url=api_host)

        if login_resp is None:
            from swanlab.sdk.internal.pkg.client.bootstrap import login_by_api_key

            login_resp = login_by_api_key(base_url=api_host + "/api", api_key=api_key)
        user_info = login_resp.get("userInfo", {}) if login_resp else {}
        username = user_info.get("username", "")
        name = user_info.get("name", "") or ""
        ctx = ApiClientContext(client=_client, web_host=web_host, api_host=api_host, username=username, name=name)
        super().__init__(ctx)

    def json(self) -> dict:
        """Api 非数据实体，返回空字典。"""
        return {}

    @staticmethod
    def _resolve_credentials(
        api_key: Optional[str],
        host: Optional[str],
        web_host: Optional[str],
    ) -> tuple[str, str, str]:
        """
        按优先级解析凭证：显式参数 > scope 登录态 > Settings（含 .netrc / 环境变量）。
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
    #  实体工厂方法
    #  - 单实体（workspace/project/run）：构造后调用 _fetch() 立即加载并返回 ok/not-ok
    #  - 列表迭代器（workspaces/projects/runs）：惰性构造，迭代时按需分页请求
    # ------------------------------------------------------------------

    def workspace(self, username: Optional[str] = None) -> Workspace:
        """
        获取工作空间信息，默认为当前登录用户的工作空间。

        :param username: 指定工作空间用户名，为 None 时使用当前登录用户
        """
        if username is None:
            username = self._ctx.username
        return Workspace(self._ctx, username=username)

    def workspaces(self, username: Optional[str] = None) -> Workspaces:
        """
        获取工作空间列表迭代器。

        :param username: 指定用户名，为 None 时使用当前登录用户
        """
        if username is None:
            username = self._ctx.username
        return Workspaces(self._ctx, username=username)

    def project(self, path: str) -> Project:
        """
        获取项目信息。

        :param path: 项目路径，格式为 'username/project-name'
        """
        return Project(self._ctx, path=path)

    def projects(
        self,
        path: str,
        sort: Optional[str] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
        page: int = 1,
        size: int = 20,
        all: bool = False,
    ) -> Projects:
        """
        获取工作空间下的项目列表迭代器。

        :param path: 工作空间名称 'username'
        :param sort: 排序方式
        :param search: 搜索关键词
        :param detail: 是否返回详细信息
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param all: 是否获取全部数据，默认 False
        """
        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        return Projects(self._ctx, path=path, query=query, detail=detail)

    def run(self, path: str) -> Experiment:
        """
        获取单个实验。

        :param path: 实验路径，格式为 'username/project/run_id'
        """

        return Experiment(self._ctx, path=path)

    def runs(
        self,
        path: str,
        filters: Optional[dict] = None,
        page: int = 1,
        size: int = 20,
        all: bool = False,
    ) -> Experiments:
        """
        获取项目下的实验列表迭代器。

        :param path: 项目路径，格式为 'username/project'
        :param filters: 筛选条件
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param all: 是否获取全部数据，默认 False
        """
        query = PaginatedQuery(page=page, size=size, all=all)
        return Experiments(self._ctx, proj_path=path, filters=filters, query=query)

    def user(self) -> User:
        return User(self._ctx)


__all__ = ["Api"]
