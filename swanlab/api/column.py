"""
@author: caddiesnew
@file: column.py
@time: 2026/4/20
@description: Column 实体类 — 实验列的查询与操作
"""

from typing import Any, Dict, Iterator, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.column import ApiColumnCsvExportType, ApiColumnType
from swanlab.api.typings.common import ApiResponseType, PaginatedQuery
from swanlab.api.utils import get_properties, parse_column_data_type, resovle_run_path, validate_column_params


class Column(BaseEntity):
    """
    表示一个 SwanLab 实验列。

    支持双模式：构造时传入 data（列表迭代注入），或 data=None（按需懒加载）。
    注意：列不支持单个获取 API，只能通过列表接口获取。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        key: str,
        column_class: Optional[str] = "CUSTOM",
        column_type: Optional[str] = None,
        data: Optional[ApiColumnType] = None,
    ) -> None:
        super().__init__(ctx)
        self._proj_path, self._run_id = resovle_run_path(path=path)
        self._key = key
        self._column_class = column_class
        self._column_type = column_type
        self._data = data
        self._project_id = None

    def _ensure_data(self) -> ApiColumnType:
        if self._data is None:
            validate_column_params(column_class=self._column_class)
            extra: Dict[str, Any] = {"search": self._key}
            if self._column_class:
                extra["class"] = self._column_class
            resp = self._get(
                f"/experiment/{self._run_id}/column",
                params={"page": 1, "size": 10, **extra},
            )
            if resp.data:
                items = resp.data.get("list", []) if isinstance(resp.data, dict) else []
                if items:
                    self._data = cast(ApiColumnType, items[0])
            if self._data is None:
                self._data = cast(ApiColumnType, {})
            self._data["run_id"] = self._run_id
        if self._project_id is None:
            resp = self._get(f"/project/{self._proj_path}")
            proj_data = resp.data if resp.ok else {}
            self._project_id = proj_data.get("cuid", "")
            self._data["project_id"] = self._project_id
        return self._data

    @property
    def project_id(self) -> str:
        if self._project_id:
            return self._project_id
        return self._ensure_data().get("project_id", "")

    @property
    def run_id(self) -> str:
        if self._run_id:
            return self._run_id
        return self._ensure_data().get("run_id", "")

    @property
    def key(self) -> str:
        res_key = self._ensure_data().get("key", "")
        if res_key and res_key != self._key:
            self._key = res_key
        return res_key

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

    def export_csv(self) -> ApiResponseType:
        """
        导出列数据为 CSV。

        :return: ApiResponseType，成功时 data 包含临时下载 URL
        """
        resp = self._get(f"/experiment/{self._run_id}/column/csv", params={"key": self.key})
        if not resp.ok:
            return resp

        data = resp.data
        if isinstance(data, list) and data:
            url = data[0].get("url", "")
            return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))
        elif isinstance(data, dict):
            url = data.get("url", "")
            return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))

        return ApiResponseType(ok=False, errmsg="Invalid response format", data=None)

    def metric(self):
        from swanlab.api.metric import Metric

        metric_type = parse_column_data_type(self.column_type)

        cur_metric = Metric(ctx=self._ctx, project_id=self.project_id, run_id=self.run_id, metric_type=metric_type)

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
        column_class: Optional[str] = None,
        column_type: Optional[str] = None,
    ) -> None:
        super().__init__(ctx)
        self._run_path = path
        self._proj_path, self._run_id = resovle_run_path(path=path)
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

    def __iter__(self) -> Iterator[Column]:
        """迭代分页获取列。"""
        extra: Dict[str, Any] = {}
        if self._column_type:
            extra["type"] = self._column_type
        if self._column_class:
            extra["class"] = self._column_class

        for item in self._paginate(
            f"/experiment/{self._run_id}/column",
            self._query,
            page_info=self._page_info,
            extra=extra,
        ):
            yield Column(
                self._ctx,
                path=self._run_path,
                key=item.get("key", ""),
                data=cast(ApiColumnType, item),
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
