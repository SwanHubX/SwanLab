"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/20
@description: SwanLab 公共查询 API 入口，面向用户的 OOP 查询接口
"""

import warnings
from typing import Any, Dict, List, Optional

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal.pkg import nrc, scope
from swanlab.sdk.internal.pkg.client import Client
from swanlab.sdk.internal.settings import settings as global_settings

from .base import ApiClientContext, BaseEntity
from .column import Column, Columns
from .experiment import Experiment, Experiments
from .project import Project, Projects
from .self_hosted import SelfHosted
from .series import Series
from .typings.common import (
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricKeyClassLiteral,
    ApiMetricKeyTypeLiteral,
    ApiVisibilityLiteral,
    PaginatedQuery,
)
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
    ) -> None:
        """
        初始化 Api 实例。

        认证优先级：
        1. 显式参数 (api_key / host)
        2. scope 登录态（进程内已调用 swanlab.login 时可用）
        3. Settings（含 .netrc / 环境变量）

        始终创建独立的 Client 实例，与 SDK 运行时单例互不干扰。

        :param api_key: API 密钥，为 None 时从 Settings / .netrc / 环境变量读取
        :param host: API 主机地址，为 None 时从 Settings 读取
        """
        # 优先从 scope 获取已有登录态（如进程内已调用 swanlab.login），直接复用凭证
        login_resp = scope.get_context("login_resp")
        api_key, api_host, web_host = self._resolve_credentials(api_key, host)
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

    @property
    def username(self) -> str:
        """当前认证用户的 username。"""
        return self._ctx.username

    @staticmethod
    def _resolve_credentials(
        api_key: Optional[str],
        host: Optional[str],
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

        if host is not None:
            api_host: str = nrc.fmt(host)
            web_host: str = api_host
        else:
            api_host = global_settings.api_host
            web_host = global_settings.web_host

        return api_key, api_host, web_host

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
        size: int = 100,
        all: bool = False,
    ) -> Projects:
        """
        获取工作空间下的项目列表迭代器。

        :param path: 工作空间名称 'username'
        :param sort: 排序方式
        :param search: 搜索关键词
        :param detail: 是否返回详细信息
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 100
        :param all: 是否获取全部数据，默认 False
        """
        validate_api_path(path, segments=1, label="workspace")
        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        return Projects(self._ctx, path=path, query=query, detail=detail)

    def create_project(
        self,
        username: Optional[str] = None,
        name: str = "",
        *,
        visibility: ApiVisibilityLiteral = "PRIVATE",
        description: Optional[str] = None,
    ) -> Optional[Project]:
        """
        在指定工作空间下创建项目。

        :param username: 工作空间用户名，为 None 时使用当前登录用户
        :param name: 项目名称 (1-100 字符，仅支持 0-9a-zA-Z-_.+)
        :param visibility: 可见性，PUBLIC 或 PRIVATE，默认 PRIVATE
        :param description: 项目描述
        """
        if username is None:
            username = self._ctx.username
        validate_api_path(username, segments=1, label="workspace")
        ws = Workspace(self._ctx, username=username)
        return ws.create_project(name, visibility=visibility, description=description)

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
    ) -> Experiments:
        """
        通过条件过滤获取项目下的实验列表。

        :param path: 项目路径，格式为 'username/project'
        :param filters: 过滤规则列表，每项为 {key, type, op, value}
        """
        validate_api_path(path, segments=2, label="project")
        if not isinstance(filters, list):
            if filters is not None:
                warnings.warn("filters must be a list, got %s — ignoring." % type(filters).__name__, stacklevel=2)
            filters = None
        return Experiments(self._ctx, path=path, filters=filters, mode="post")

    def runs_get(
        self,
        path: str,
        page: int = 1,
        size: int = 100,
        all: bool = False,
    ) -> Experiments:
        """
        通过分页获取项目下的实验列表。

        :param path: 项目路径，格式为 'username/project'
        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 100
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
        size: int = 100,
        search: Optional[str] = None,
        column_class: ApiColumnClassLiteral = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        all: bool = False,
    ) -> Columns:
        """
        .. deprecated::
            Use :meth:`series` instead. Column metadata (name/class/type) is no longer
            maintained by the House backend; ``series()`` provides cursor-paginated
            key listing directly from House.

        List columns under an experiment (paginated, with optional fuzzy search).

        The ``search`` parameter performs **fuzzy matching** (case-insensitive ``contains``)
        on the column ``name`` field.

        :param path: Experiment path, format: ``'username/project/run_id'``
        :param page: Page number, default 1
        :param size: Page size, default 100
        :param search: Fuzzy search keyword (matches column **name**, not key)
        :param column_class: Column class, ``CUSTOM`` or ``SYSTEM``, default ``CUSTOM``
        :param column_type: Column data type, e.g. ``FLOAT``, ``STRING``, ``IMAGE``
        :param all: If True, fetch all pages, default False
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
        .. deprecated::
            Use :meth:`series` instead. Column metadata (name/class/type) is no longer
            maintained by the House backend; use ``series()`` for per-key access
            (e.g. ``series(keys=["loss"])``).

        Get a single column by key (fuzzy search, first match).

        Performs fuzzy search (``contains`` on ``name``) and returns the first matching
        column. If multiple columns share a similar name, the first one (ordered by
        ``id DESC``) is returned.

        :param path: Experiment path, format: ``'username/project/run_id'``
        :param key: Column key to search, e.g. ``"loss"``, ``"acc"``
        :param column_class: Column class, ``CUSTOM`` or ``SYSTEM``, default ``CUSTOM``
        :param column_type: Column data type, e.g. ``FLOAT``, ``STRING``, ``IMAGE``
        """
        validate_api_path(path, segments=3, label="run")
        validate_non_empty_string(key, label="column key")
        return Column(self._ctx, path=path, key=key, column_class=column_class, column_type=column_type)

    def series(
        self,
        path: str,
        metric_type: ApiMetricKeyTypeLiteral = "SCALAR",
        limit: int = 2000,
        all: bool = False,
        metric_class: ApiMetricKeyClassLiteral = "CUSTOM",
    ) -> Series:
        """
        Cursor-paginated listing of metric keys (preferred over deprecated ``columns()``).

        :param path: Experiment path, format: ``'username/project/run_id'``
        :param metric_type: ``SCALAR`` (default) or ``MEDIA``
        :param limit: Max keys per page (1..2000)
        :param all: If False (default), return only one page (up to ``limit`` keys) with ``hasMore``
            indicating more keys exist. If True, iterate all pages.
        :param metric_class: ``CUSTOM`` (default) or ``SYSTEM``. Post-filter keys by class.
            SCALAR keys with ``__swanlab__`` prefix are SYSTEM; all MEDIA keys are CUSTOM.
        """
        validate_api_path(path, segments=3, label="run")
        exp = Experiment(self._ctx, path=path)
        return exp.series(metric_type=metric_type, limit=limit, all=all, metric_class=metric_class)

    # -------
    # 私有化相关接口
    # --------
    def self_hosted(self) -> SelfHosted:
        return SelfHosted(self._ctx)


__all__ = ["Api", "Series"]
