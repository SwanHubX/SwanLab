"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/20
@description: SwanLab public query API entry — OOP query interface for users
"""

import warnings
from typing import Any, Dict, List, Optional

from typing_extensions import deprecated

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
    SwanLab public query API entry.

    Usage::

        from swanlab import Api

        api = Api()                              # auto-read credentials from .netrc
        api = Api(api_key="...", host="...")      # explicit credentials

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
        Initialize an Api instance.

        Credential resolution order:
        1. Explicit parameters (``api_key`` / ``host``)
        2. In-process login state (available when ``swanlab.login`` has been called)
        3. Settings (including ``.netrc`` / environment variables)

        :param api_key: API key; if None, resolved per the order above
        :param host: API host address; if None, uses default configuration
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
        """Return an empty dict."""
        return {}

    @property
    def username(self) -> str:
        """The authenticated user's username."""
        return self._ctx.username

    @staticmethod
    def _resolve_credentials(
        api_key: Optional[str],
        host: Optional[str],
    ) -> tuple[str, str, str]:
        """
        Resolve credentials by priority: explicit params > in-process login state > Settings (.netrc / env vars).
        Returns (api_key, api_host, web_host).
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
        Get workspace info, defaulting to the current logged-in user.

        :param username: Workspace username; None uses the current logged-in user
        """
        if username is None:
            username = self._ctx.username
        validate_api_path(username, segments=1, label="workspace")
        return Workspace(self._ctx, username=username)

    def workspaces(self, username: Optional[str] = None) -> Workspaces:
        """
        Get a workspace list iterator.

        :param username: Username; None uses the current logged-in user
        """
        if username is None:
            username = self._ctx.username
        validate_api_path(username, segments=1, label="workspace")
        return Workspaces(self._ctx, username=username)

    def project(self, path: str) -> Project:
        """
        Get project info.

        :param path: Project path, format: ``'username/project-name'``
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
        Get a project list iterator under a workspace.

        :param path: Workspace name, e.g. ``'username'``
        :param sort: Sort field
        :param search: Search keyword
        :param detail: Whether to return detailed info
        :param page: Page number, default 1
        :param size: Page size, default 100
        :param all: If True, fetch all pages, default False
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
        Create a project under the specified workspace.

        :param username: Workspace username; None uses the current logged-in user
        :param name: Project name (1-100 chars, only 0-9a-zA-Z-_.+ supported)
        :param visibility: Visibility, ``PUBLIC`` or ``PRIVATE``, default ``PRIVATE``
        :param description: Project description
        """
        if username is None:
            username = self._ctx.username
        validate_api_path(username, segments=1, label="workspace")
        ws = Workspace(self._ctx, username=username)
        return ws.create_project(name, visibility=visibility, description=description)

    def run(self, path: str) -> Experiment:
        """
        Get a single experiment.

        :param path: Experiment path, format: ``'username/project/run_id'``
        """
        validate_api_path(path, segments=3, label="run")
        return Experiment(self._ctx, path=path)

    def runs(
        self,
        path: str,
        filters: Optional[List[Dict[str, Any]]] = None,
    ) -> Experiments:
        """
        Get a filtered experiment list under a project.

        :param path: Project path, format: ``'username/project'``
        :param filters: Filter rules list, each item is ``{key, type, op, value}``
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
        Get a paginated experiment list under a project.

        :param path: Project path, format: ``'username/project'``
        :param page: Page number, default 1
        :param size: Page size, default 100
        :param all: If True, fetch all pages, default False
        """
        validate_api_path(path, segments=2, label="project")
        query = PaginatedQuery(page=page, size=size, all=all)
        return Experiments(self._ctx, path=path, query=query, mode="get")

    def user(self) -> User:
        return User(self._ctx)

    @deprecated("Use `series()` method instead.")
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
        List columns under an experiment (paginated, with optional fuzzy search).

        :param path: Experiment path, format: ``'username/project/run_id'``
        :param page: Page number, default 1
        :param size: Page size, default 100
        :param search: Fuzzy search keyword
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

    @deprecated("Use `series()` method instead.")
    def column(
        self,
        path: str,
        key: str,
        column_class: Optional[ApiColumnClassLiteral] = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
    ) -> Column:
        """
        Get a single column by key (fuzzy search, first match).

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
        metric_class: ApiMetricKeyClassLiteral = "CUSTOM",
        search: str = "",
    ) -> Series:
        """
        List all metric keys for the given experiment.

        :param path: Experiment path, format: ``'username/project/run_id'``
        :param metric_type: ``"SCALAR"`` (default) or ``"MEDIA"``
        :param metric_class: ``"CUSTOM"`` (default) or ``"SYSTEM"`` — filter keys by class.
        :param search: Fuzzy search filter — case-insensitive substring match on key names.
        :returns: :class:`Series`
        """
        validate_api_path(path, segments=3, label="run")
        exp = Experiment(self._ctx, path=path)
        return exp.series(metric_type=metric_type, metric_class=metric_class, search=search)

    # -------
    # 私有化相关接口
    # --------
    def self_hosted(self) -> SelfHosted:
        return SelfHosted(self._ctx)


__all__ = ["Api", "Series"]
