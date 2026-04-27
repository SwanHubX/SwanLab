"""
@author: caddiesnew
@file: column.py
@time: 2026/4/20
@description: Column 实体类 — 实验列的查询与操作
"""

from typing import Any, Callable, Dict, Iterator, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.column import ApiColumnType
from swanlab.api.typings.common import (
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricColumnTypeLiteral,
    ApiResponseType,
    PaginatedQuery,
)
from swanlab.api.utils import get_properties, parse_column_data_type, resolve_run_path, validate_column_params


class Column(BaseEntity):
    """
    表示一个 SwanLab 实验列。

    支持双模式：构造时传入 data（列表迭代注入），或 data=None（按需懒加载）。
    注意：列不支持单个获取 API，只能通过列表接口获取。
    """

    @staticmethod
    def _resolve_cuid(entity: BaseEntity, path: str, fallback: str = "") -> str:
        resp = entity._get(path)
        data = resp.data if resp.ok and isinstance(resp.data, dict) else {}
        return data.get("cuid", "") or fallback

    @staticmethod
    def _resolve_run_cuid(entity: BaseEntity, project_path: str, run_slug: str) -> str:
        return Column._resolve_cuid(entity, f"/project/{project_path}/runs/{run_slug}", fallback=run_slug)

    @staticmethod
    def _resolve_project_cuid(entity: BaseEntity, project_path: str) -> str:
        return Column._resolve_cuid(entity, f"/project/{project_path}")

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        key: str,
        column_class: Optional[ApiColumnClassLiteral] = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        data: Optional[ApiColumnType] = None,
        project_id: Optional[str] = None,
        run_id: Optional[str] = None,
        project_id_getter: Optional[Callable[[], str]] = None,
    ) -> None:
        super().__init__(ctx)
        self._proj_path, self._run_slug = resolve_run_path(path=path)
        self._key = key
        self._column_class = column_class
        self._column_type = column_type
        self._data = data
        self._project_id = project_id or (data or {}).get("project_id", "") or None
        self._run_id = run_id or (data or {}).get("run_id", "") or ""
        self._project_id_getter = project_id_getter

    def _ensure_data(self) -> ApiColumnType:
        if self._data is None:
            validate_column_params(column_class=self._column_class, column_type=self._column_type)
            extra: Dict[str, Any] = {"search": self._key}
            run_id = self._ensure_run_id()
            if self._column_class:
                extra["class"] = self._column_class
            resp = self._get(
                f"/experiment/{run_id}/column",
                params={"page": 1, "size": 10, **extra},
            )
            if resp.data:
                items = resp.data.get("list", []) if isinstance(resp.data, dict) else []
                if items:
                    self._data = cast(ApiColumnType, items[0])
            if self._data is None:
                self._data = cast(ApiColumnType, {})
        self._data.setdefault("run_id", self._ensure_run_id())
        return self._data

    def _ensure_run_id(self) -> str:
        if self._run_id:
            return self._run_id
        if self._data and (run_id := self._data.get("run_id", "")):
            self._run_id = run_id
            return run_id
        self._run_id = Column._resolve_run_cuid(self, self._proj_path, self._run_slug)
        if self._data is not None:
            self._data["run_id"] = self._run_id
        return self._run_id

    def _ensure_project_id(self) -> str:
        if self._project_id:
            return self._project_id
        if self._data and (project_id := self._data.get("project_id", "")):
            self._project_id = project_id
            return project_id
        if self._project_id_getter is not None:
            self._project_id = self._project_id_getter()
            if self._data is not None:
                self._data["project_id"] = self._project_id
            return self._project_id
        if self._project_id is None:
            self._project_id = Column._resolve_project_cuid(self, self._proj_path)
            if self._data is not None:
                self._data["project_id"] = self._project_id
        return self._project_id

    @property
    def project_id(self) -> str:
        return self._ensure_project_id()

    @property
    def run_id(self) -> str:
        return self._ensure_run_id()

    @property
    def key(self) -> str:
        if self._key:
            return self._key
        return self._ensure_data().get("key", "")

    @property
    def name(self) -> str:
        """列的显示名称，默认为 key 的值。"""
        return self._ensure_data().get("name", "")

    @property
    def column_class(self) -> str:
        """列的分类：CUSTOM 或 SYSTEM。"""
        return self._ensure_data().get("class", "")

    @property
    def column_type(self) -> str:
        """列的数据类型，如 FLOAT、STRING、IMAGE 等。"""
        return self._ensure_data().get("type", "")

    @property
    def created_at(self) -> int:
        """列的创建时间戳。"""
        return self._ensure_data().get("createdAt", 0)

    @property
    def error(self) -> Optional[Dict[str, Any]]:
        """列的错误信息。"""
        return self._ensure_data().get("error", {})

    def metric(
        self,
        sample: int = 1500,
        metric_type: ApiMetricColumnTypeLiteral = "SCALAR",
        ignore_timestamp: bool = False,
        media_step: Optional[int] = None,
    ) -> Dict[str, Any]:
        from swanlab.api.metric import Metric

        metric_type = parse_column_data_type(self.column_type)
        metric = Metric(
            ctx=self._ctx,
            project_id=self.project_id,
            run_id=self.run_id,
            key=self.key,
            sample=sample,
            metric_type=metric_type,
            ignore_timestamp=ignore_timestamp,
            media_step=media_step,
        )
        return metric.json()

    def export_csv(self) -> ApiResponseType:
        from swanlab.api.metric import Metric

        metric_type = parse_column_data_type(self.column_type)
        if metric_type != "SCALAR":
            err_msg = "export_csv() only support SCALAR metric_type"
            return ApiResponseType(ok=False, errmsg=err_msg, data=None)
        metric = Metric(
            ctx=self._ctx,
            project_id=self.project_id,
            run_id=self.run_id,
            key=self.key,
            metric_type=metric_type,
        )
        return metric.export_csv()

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Columns(BaseEntity):
    """
    实验下列集合的分页迭代器。

    用法::

        # 获取所有列
        for column in experiment.columns():
            print(column.name, column.data_type)

        # 分页获取列（支持搜索）
        for column in experiment.columns(page=1, size=20, search="loss"):
            print(column.name)

        # 获取全部列（自动翻页）
        for column in experiment.columns(all=True):
            print(column.name)
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        query: PaginatedQuery,
        column_class: Optional[ApiColumnClassLiteral] = None,
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        project_id: Optional[str] = None,
        run_id: Optional[str] = None,
        project_id_getter: Optional[Callable[[], str]] = None,
    ) -> None:
        super().__init__(ctx)
        self._run_path = path
        self._proj_path, self._run_slug = resolve_run_path(path=path)
        self._run_id = run_id or ""
        self._project_id = project_id
        self._project_id_getter = project_id_getter
        self._query = query
        # 校验 column_type 和 column_class 的合法性
        validate_column_params(column_type=column_type, column_class=column_class)
        self._column_class = column_class
        self._column_type = column_type
        self._page_info: Dict[str, Any] = {
            "page": query.page,
            "size": query.size,
            "total": 0,
            "pages": 0,
            "list": [],
        }

    def _ensure_run_id(self) -> str:
        if self._run_id:
            return self._run_id
        self._run_id = Column._resolve_run_cuid(self, self._proj_path, self._run_slug)
        return self._run_id

    def _ensure_project_id(self) -> str:
        if self._project_id:
            return self._project_id
        if self._project_id_getter is not None:
            self._project_id = self._project_id_getter()
            return self._project_id
        self._project_id = Column._resolve_project_cuid(self, self._proj_path)
        return self._project_id

    def __iter__(self) -> Iterator[Column]:
        """迭代分页获取列。"""
        extra: Dict[str, Any] = {}
        run_id = self._ensure_run_id()
        if self._column_type:
            extra["type"] = self._column_type
        if self._column_class:
            extra["class"] = self._column_class

        for item in self._paginate(
            f"/experiment/{run_id}/column",
            self._query,
            page_info=self._page_info,
            extra=extra,
        ):
            data = {**item, "run_id": run_id}
            yield Column(
                self._ctx,
                path=self._run_path,
                key=item.get("key", ""),
                data=cast(ApiColumnType, data),
                run_id=run_id,
                project_id_getter=self._ensure_project_id,
            )

    @property
    def total(self) -> int:
        """获取总数（触发一次请求）。"""
        # 触发一次迭代来获取总数
        if self._page_info["total"] == 0:
            try:
                next(iter(self))
            except StopIteration:
                pass
        return self._page_info["total"]

    def json(self) -> Dict[str, Any]:
        self._page_info["list"] = [c.json() for c in self]
        return self._page_info
