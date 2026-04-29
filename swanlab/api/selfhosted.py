"""
@author: caddiesnew
@file: selfhosted.py
@time: 2026/4/20
@description: SelfHosted 实体类 — 私有化部署实例的查询与管理
"""

from typing import Any, Dict, Iterator, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import ApiResponseType, PaginatedQuery
from swanlab.api.typings.selfhosted import ApiLicensePlanLiteral, ApiSelfHostedInfoType
from swanlab.api.utils import get_properties, validate_non_empty_string


class SelfHosted(BaseEntity):
    """
    表示一个 SwanLab 私有化部署实例。

    支持双模式：构造时传入 data，或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        data: Optional[ApiSelfHostedInfoType] = None,
    ) -> None:
        super().__init__(ctx)
        self._data = data

    def _ensure_data(self) -> ApiSelfHostedInfoType:
        if self._data is None:
            resp = self._get("/self_hosted/info")
            self._data = resp.data if resp.ok and resp.data else cast(ApiSelfHostedInfoType, {})
        return self._data

    @property
    def enabled(self) -> bool:
        return self._ensure_data().get("enabled", False)

    @property
    def expired(self) -> bool:
        return self._ensure_data().get("expired", False)

    @property
    def root(self) -> bool:
        return self._ensure_data().get("root", False)

    @property
    def plan(self) -> ApiLicensePlanLiteral:
        return self._ensure_data().get("plan", "free")

    @property
    def seats(self) -> int:
        return self._ensure_data().get("seats", 0)

    # ================================
    #   权限校验
    # ================================

    @staticmethod
    def validate_expire(info: ApiSelfHostedInfoType) -> None:
        if info.get("expired", True):
            raise ValueError("SwanLab self-hosted instance has expired.")

    @staticmethod
    def validate_root(info: ApiSelfHostedInfoType) -> None:
        SelfHosted.validate_expire(info)
        if not info.get("root", False):
            raise ValueError("You don't have permission to perform this action. Please login as a root user.")

    # ================================
    #   管理操作（root 限定）
    # ================================

    def create_user(self, username: str, password: str) -> ApiResponseType:
        """
        添加用户（私有化管理员限定）。

        :param username: 待创建用户名
        :param password: 待创建用户密码
        """
        SelfHosted.validate_root(self._ensure_data())
        validate_non_empty_string(username, label="username")
        validate_non_empty_string(password, label="password")
        data = {"users": [{"username": username, "password": password}]}
        return self._post("/self_hosted/users", data=data)

    def get_users(self, page: int = 1, size: int = 20, all: bool = False) -> Iterator[dict]:
        """
        分页获取用户（管理员限定）。

        :param page: 起始页码，默认 1
        :param size: 每页大小，默认 20
        :param all: 是否获取全部数据，默认 False
        """
        SelfHosted.validate_root(self._ensure_data())
        query = PaginatedQuery(page=page, size=size, all=all)
        page_info: Dict[str, Any] = {"total": 0, "pages": 0}
        yield from self._paginate("/self_hosted/users", query, page_info=page_info)

    def get_projects(
        self,
        page: int = 1,
        size: int = 20,
        search: Optional[str] = None,
        sort: Optional[str] = None,
        state: Optional[str] = None,
        creator: Optional[str] = None,
        group: Optional[str] = None,
        all: bool = False,
    ) -> Iterator[dict]:
        """
        分页获取所有项目（管理员限定）。

        :param page: 起始页码，默认 1
        :param size: 每页大小，默认 20
        :param search: 搜索关键词
        :param sort: 排序字段，update(默认) / create / name
        :param state: 实验状态过滤，RUNNING / FINISHED
        :param creator: 创建者 username
        :param group: 组织空间 username
        :param all: 是否获取全部数据，默认 False
        """
        SelfHosted.validate_root(self._ensure_data())
        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        page_info: Dict[str, Any] = {"total": 0, "pages": 0}
        yield from self._paginate(
            "/self_hosted/projects",
            query,
            page_info=page_info,
            extra={"state": state, "creator": creator, "group": group},
        )

    def get_groups(
        self,
        page: int = 1,
        size: int = 20,
        search: Optional[str] = None,
        type: Optional[str] = None,
        sort: Optional[str] = None,
        all: bool = False,
    ) -> Iterator[dict]:
        """
        分页获取所有空间（管理员限定）。

        :param page: 起始页码，默认 1
        :param size: 每页大小，默认 20
        :param search: 搜索关键词
        :param type: 空间类型过滤，PERSON / TEAM
        :param sort: 排序字段，update(默认) / create / name
        :param all: 是否获取全部数据，默认 False
        """
        SelfHosted.validate_root(self._ensure_data())
        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        page_info: Dict[str, Any] = {"total": 0, "pages": 0}
        yield from self._paginate(
            "/self_hosted/groups",
            query,
            page_info=page_info,
            extra={"type": type},
        )

    def get_usage_summary(self) -> ApiResponseType:
        """获取系统内汇总信息（root 限定）"""
        SelfHosted.validate_root(self._ensure_data())
        return self._get("/self_hosted/summary")

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
