"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/20
@description: SwanLab 公共查询 API 入口，面向用户的 OOP 查询接口
"""

from typing import Any, Dict, List, Optional

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.pkg import nrc, scope
from swanlab.sdk.internal.pkg.client import Client
from swanlab.sdk.internal.settings import settings as global_settings

from .base import ApiClientContext, BaseEntity
from .column import Column, Columns
from .experiment import Experiment, Experiments
from .project import Project, Projects
from .selfhosted import SelfHosted
from .typings.common import ApiColumnClassLiteral, ApiColumnDataTypeLiteral, PaginatedQuery
from .user import User
from .utils import validate_api_path, validate_non_empty_string
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
        if not isinstance(api_key, str) or not api_key.strip():
            raise AuthenticationError("No API key found. Please login with `swanlab login` or pass api_key parameter.")
        api_key = api_key.strip()

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
        validate_api_path(username, segments=1, label="workspace")
        return Workspace(self._ctx, username=username)

    def workspaces(self, username: Optional[str] = None) -> Workspaces:
        """
        获取工作空间列表迭代器。

        :param username: 指定用户名，为 None 时使用当前登录用户
        """
        if username is None:
            username = self._ctx.username
        validate_api_path(username, segments=1, label="workspace")
        return Workspaces(self._ctx, username=username)

    def project(self, path: str) -> Project:
        """
        获取项目信息。

        :param path: 项目路径，格式为 'username/project-name'
        """
        validate_api_path(path, segments=2, label="project")
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
        validate_api_path(path, segments=1, label="workspace")
        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        return Projects(self._ctx, path=path, query=query, detail=detail)

    def run(self, path: str) -> Experiment:
        """
        获取单个实验。

        :param path: 实验路径，格式为 'username/project/run_id'
        """
        validate_api_path(path, segments=3, label="run")
        return Experiment(self._ctx, path=path)

    def runs(
        self,
        path: str,
        filters: Optional[List[Dict[str, Any]]] = None,
        groups: Optional[List[Dict[str, Any]]] = None,
        sorts: Optional[List[Dict[str, Any]]] = None,
    ) -> Experiments:
        """
        通过条件过滤获取项目下的实验列表。

        :param path: 项目路径，格式为 'username/project'
        :param filters: 过滤规则列表，每项为 {key, type, op, value}
        :param groups: 分组规则列表，每项为 {key, type}
        :param sorts: 排序规则列表，每项为 {key, type, order}
        """
        validate_api_path(path, segments=2, label="project")
        return Experiments(self._ctx, path=path, filters=filters, groups=groups, sorts=sorts, mode="post")

    def runs_get(
        self,
        path: str,
        page: int = 1,
        size: int = 20,
        all: bool = False,
    ) -> Experiments:
        """
        通过分页获取项目下的实验列表。

        :param path: 项目路径，格式为 'username/project'
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param all: 是否获取全部数据，默认 False
        """
        validate_api_path(path, segments=2, label="project")
        query = PaginatedQuery(page=page, size=size, all=all)
        return Experiments(self._ctx, path=path, query=query, mode="get")

    def user(self) -> User:
        return User(self._ctx)

    def columns(
        self,
        path: str,
        page: int = 1,
        size: int = 20,
        search: Optional[str] = None,
        column_class: ApiColumnClassLiteral = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        all: bool = False,
    ) -> Columns:
        """
        获取实验下的列列表（分页查询，支持搜索）。

        :param path: 实验路径，格式为 'username/project/run_id'
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param search: 搜索关键词，搜索的是列的 name
        :param column_class: 列的分类，CUSTOM 或 SYSTEM, 默认为 CUSTOM
        :param column_type: 列的类型，如 FLOAT、STRING、IMAGE 等
        :param all: 是否获取全部数据，默认 False
        """
        validate_api_path(path, segments=3, label="run")
        query = PaginatedQuery(page=page, size=size, search=search, all=all)
        return Columns(
            self._ctx,
            path=path,
            query=query,
            column_type=column_type,
            column_class=column_class,
        )

    def column(
        self,
        path: str,
        key: str,
        column_class: Optional[ApiColumnClassLiteral] = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
    ) -> Column:
        """
        获取单个列（通过搜索 key 匹配）。

        :param path: 实验路径，格式为 'username/project/run_id'
        :param key: 列的键名, 输入不完整则模糊匹配 name 为首个 key.
        :param column_class: 列的分类，CUSTOM 或 SYSTEM，默认 CUSTOM
        :param column_type: 列的类型，如 FLOAT、STRING、IMAGE 等，默认为 None
        """
        validate_api_path(path, segments=3, label="run")
        validate_non_empty_string(key, label="column key")
        return Column(self._ctx, path=path, key=key, column_class=column_class, column_type=column_type)

    # -------
    # 私有化相关接口
    # --------
    def self_hosted(self) -> SelfHosted:
        return SelfHosted(self._ctx)


__all__ = ["Api"]
