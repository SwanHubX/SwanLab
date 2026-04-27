"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: Project 实体类 — 单个项目的查询与操作
"""

from typing import Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import PaginatedQuery
from swanlab.api.typings.project import ApiProjectCountType, ApiProjectLabelType, ApiProjectType
from swanlab.api.utils import get_properties


class Project(BaseEntity):
    """
    表示一个 SwanLab 项目。

    支持双模式：构造时传入 data（列表迭代注入），或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        data: Optional[ApiProjectType] = None,
    ) -> None:
        super().__init__(ctx)
        self._path = path
        self._data = data

    def _ensure_data(self) -> ApiProjectType:
        if self._data is None:
            resp = self._get(f"/project/{self._path}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiProjectType, {})
        return self._data

    @property
    def project_id(self) -> str:
        return self._ensure_data().get("cuid", "")

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def path(self) -> str:
        return self._ensure_data().get("path", "")

    @property
    def url(self) -> str:
        return self._build_web_url(f"@{self.path}")

    @property
    def description(self) -> str:
        return self._ensure_data().get("description", "")

    @property
    def visibility(self) -> str:
        return self._ensure_data().get("visibility", "PRIVATE")

    @property
    def created_at(self) -> str:
        return self._ensure_data().get("createdAt", "")

    @property
    def updated_at(self) -> str:
        return self._ensure_data().get("updatedAt", "")

    @property
    def labels(self) -> List[ApiProjectLabelType]:
        return [label for label in self._ensure_data().get("projectLabels", [])]

    @property
    def count(self) -> ApiProjectCountType:
        return cast(ApiProjectCountType, self._ensure_data().get("_count", {}))

    def runs(
        self,
        filters: Optional[List[Dict[str, Any]]] = None,
        groups: Optional[List[Dict[str, Any]]] = None,
        sorts: Optional[List[Dict[str, Any]]] = None,
    ):
        """
        获取项目下的实验列表（POST 模式，支持复杂过滤）。

        :param filters: 过滤规则列表，每项为 {key, type, op, value}
        :param groups: 分组规则列表，每项为 {key, type}
        :param sorts: 排序规则列表，每项为 {key, type, order}
        """
        from swanlab.api.experiment import Experiments

        return Experiments(self._ctx, path=self.path, filters=filters, groups=groups, sorts=sorts, mode="post")

    def runs_get(
        self,
        page: int = 1,
        size: int = 20,
        all: bool = False,
    ):
        """
        获取项目下的实验列表（GET 模式，标准分页，返回精简信息）。

        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param all: 是否获取全部数据，默认 False
        """
        from swanlab.api.experiment import Experiments

        query = PaginatedQuery(page=page, size=size, all=all)
        return Experiments(self._ctx, path=self.path, query=query, mode="get")

    def delete(self) -> bool:
        """删除此项目。"""
        resp = self._delete(f"/project/{self.path}")
        return resp.ok

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Projects(BaseEntity):
    """
    工作空间下项目集合的分页迭代器。

    用法::

        for project in api.projects("username"):
            print(project.name)
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        query: Optional[PaginatedQuery] = None,
        detail: Optional[bool] = True,
    ) -> None:
        super().__init__(ctx)
        self._path = path
        self._query = query or PaginatedQuery()
        self._detail = detail
        self._page_info: Dict[str, Any] = {
            "page": self._query.page,
            "size": self._query.size,
            "total": 0,
            "pages": 0,
            "list": [],
        }

    def __iter__(self) -> Iterator[Project]:
        for item in self._paginate(
            f"/project/{self._path}",
            self._query,
            page_info=self._page_info,
            extra={"detail": self._detail},
        ):
            yield Project(
                self._ctx,
                path=str(item.get("path", "")),
                data=cast(ApiProjectType, item),
            )

    def json(self) -> Dict[str, Any]:
        self._page_info["list"] = [p.json() for p in self]
        return self._page_info
