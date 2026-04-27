"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from typing import Any, Dict, Iterator, List, Optional, Union, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import (
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricLogLevelLiteral,
    PaginatedQuery,
)
from swanlab.api.typings.experiment import (
    ApiExperimentLabelType,
    ApiExperimentProfileType,
    ApiExperimentType,
)
from swanlab.api.typings.user import ApiUserType
from swanlab.api.utils import (
    get_properties,
    resolve_run_path,
    validate_filter,
    validate_group,
    validate_sort,
    validate_update_active,
)


class Experiment(BaseEntity):
    """
    表示一个 SwanLab 实验（完整信息，通过 POST /runs/shows 或单实验详情接口获取）。

    支持双模式：构造时传入 data，或 data=None（按需懒加载）。
    构造时从 data 中提取 _cuid 缓存，避免 _ensure_data 与 id 属性的循环调用。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        data: Optional[ApiExperimentType] = None,
    ) -> None:
        super().__init__(ctx)
        self._proj_path, self._run_slug = resolve_run_path(path=path)
        self._cuid: str = (data or {}).get("cuid", "") or self._run_slug
        self._data: Optional[ApiExperimentType] = data
        self._project_id = ""

    def _refresh_cuid(self) -> None:
        if self._data:
            self._cuid = self._data.get("cuid", "") or self._cuid

    def _ensure_data(self) -> ApiExperimentType:
        if self._data is None:
            resp = self._get(f"/project/{self._proj_path}/runs/{self._cuid}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiExperimentType, {})
        self._refresh_cuid()
        assert self._data is not None
        return self._data

    def _ensure_project_id(self) -> str:
        if self._project_id:
            return self._project_id
        if self._data and (project_id := str(self._data.get("project_id") or "")):
            self._project_id = project_id
            return project_id
        if not self._project_id:
            resp = self._get(f"/project/{self._proj_path}")
            proj_data = resp.data if resp.ok else {}
            self._project_id = str(proj_data.get("cuid") or "")
            if self._data is not None:
                self._data["project_id"] = self._project_id
        return self._project_id

    @property
    def project_id(self) -> str:
        return self._ensure_project_id()

    @property
    def run_id(self) -> str:
        self._ensure_data()
        return self._cuid

    def _run_url_ref(self) -> str:
        data = self._ensure_data()
        return data.get("slug") or self._run_slug or self.run_id

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def description(self) -> str:
        return self._ensure_data().get("description", "")

    @property
    def type(self) -> str:
        return self._ensure_data().get("type", "")

    @property
    def state(self) -> str:
        return self._ensure_data().get("state", "")

    @property
    def url(self) -> str:
        return self._build_web_url(f"@{self._proj_path}/runs/{self._run_url_ref()}/chart")

    @property
    def show(self) -> bool:
        return self._ensure_data().get("show", True)

    @property
    def labels(self) -> List[ApiExperimentLabelType]:
        return [label for label in self._ensure_data().get("labels", [])]

    @property
    def group(self) -> str:
        return self._ensure_data().get("cluster", "")

    @property
    def job_type(self) -> str:
        return self._ensure_data().get("job", "")

    @property
    def user(self) -> ApiUserType:
        user_data = self._ensure_data().get("user", {})
        return cast(ApiUserType, user_data)

    @property
    def created_at(self) -> str:
        return self._ensure_data().get("createdAt", "")

    @property
    def finished_at(self) -> str:
        return self._ensure_data().get("finishedAt", "")

    @property
    def profile(self) -> ApiExperimentProfileType:
        """Experiment profile containing config, metadata, requirements, and conda."""
        data = self._ensure_data()
        if "profile" not in data and self._cuid:
            resp = self._get(f"/project/{self._proj_path}/runs/{self._cuid}")
            if resp.ok and resp.data:
                self._data = resp.data
                self._refresh_cuid()
                data = self._data
        return cast(ApiExperimentProfileType, self._ensure_data().get("profile", {}))

    def column(
        self,
        key: str,
        column_class: Optional[ApiColumnClassLiteral] = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = "FLOAT",
    ):
        """
        获取实验下指定 key 的单个列。

        :param key: 列的 key，如 "loss"、"acc"
        :param column_class: 列的分类，CUSTOM 或 SYSTEM
        :param column_type: 列的数据类型，如 FLOAT、STRING、IMAGE 等
        """
        from swanlab.api.column import Column

        run_id = self.run_id
        return Column(
            self._ctx,
            path=f"{self._proj_path}/{run_id}",
            key=key,
            column_class=column_class,
            column_type=column_type,
            run_id=run_id,
            project_id_getter=lambda: self.project_id,
        )

    def metrics(
        self,
        keys: List[str],
        sample: int = 1500,
        ignore_timestamp: bool = True,
        all: bool = False,
    ) -> Dict[str, Any]:
        from swanlab.api.metric import Metrics

        return Metrics(
            ctx=self._ctx,
            project_id=self.project_id,
            run_id=self.run_id,
            keys=keys,
            sample=sample,
            metric_type="SCALAR",
            ignore_timestamp=ignore_timestamp,
            all=all,
        ).json()

    def medias(
        self,
        keys: List[str],
        step: Optional[int] = 0,
        all: bool = False,
    ) -> Dict[str, Any]:
        from swanlab.api.metric import Metrics

        return Metrics(
            ctx=self._ctx,
            project_id=self.project_id,
            run_id=self.run_id,
            keys=keys,
            metric_type="MEDIA",
            media_step=step,
            all=all,
        ).json()

    def logs(
        self,
        offset: Optional[int] = 0,
        level: ApiMetricLogLevelLiteral = "INFO",
        ignore_timestamp: bool = True,
    ) -> Dict[str, Any]:
        from swanlab.api.metric import Metric

        logs = Metric(
            ctx=self._ctx,
            project_id=self.project_id,
            run_id=self.run_id,
            key="LOG",
            log_offset=offset,
            log_level=level,
            metric_type="LOG",
            ignore_timestamp=ignore_timestamp,
        )
        return logs.json()

    def columns(
        self,
        page: int = 1,
        size: int = 20,
        search: Optional[str] = None,
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        column_class: Optional[ApiColumnClassLiteral] = None,
        all: bool = False,
    ):
        """
        获取实验下的列列表（分页查询，支持搜索）。

        :param page: 起始页码，默认 1
        :param size: 每页数量，默认 20
        :param search: 搜索关键词，搜索的是列的 name
        :param column_type: 列的类型，如 FLOAT、STRING、IMAGE 等
        :param column_class: 列的分类，CUSTOM 或 SYSTEM
        :param all: 是否获取全部数据，默认 False
        """
        from swanlab.api.column import Columns

        query = PaginatedQuery(page=page, size=size, search=search, all=all)
        run_id = self.run_id
        return Columns(
            self._ctx,
            path=f"{self._proj_path}/{run_id}",
            query=query,
            column_type=column_type,
            column_class=column_class,
            run_id=run_id,
            project_id_getter=lambda: self.project_id,
        )

    def delete(self) -> bool:
        """删除此实验。"""
        resp = self._delete(f"/project/{self._proj_path}/runs/{self.run_id}")
        return resp.ok

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


def _flatten_runs(runs: Union[list, Dict]) -> list:
    """展开分组后的实验数据，返回一个包含所有实验的列表。"""
    if isinstance(runs, dict):
        return [item for v in runs.values() for item in _flatten_runs(v)]
    if isinstance(runs, list):
        return list(runs)
    return [runs]


class Experiments(BaseEntity):
    """
    项目下实验集合的迭代器。

    支持两种模式：
    - POST 模式（默认）：通过 /runs/shows 接口获取，支持复杂过滤，不支持分页
    - GET 模式：通过 /runs 接口获取，支持标准分页，返回精简信息

    用法::

        # POST 复杂过滤
        for run in api.runs(path="username/project"):
            print(run.name)

        # GET 分页
        for run in api.list_runs_simple(path="username/project"):
            print(run.name, run.state)
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        path: str,
        filters: Optional[List[Dict[str, Any]]] = None,
        groups: Optional[List[Dict[str, Any]]] = None,
        sorts: Optional[List[Dict[str, Any]]] = None,
        query: Optional[PaginatedQuery] = None,
        mode: str = "post",
    ) -> None:
        super().__init__(ctx)
        self._proj_path = path
        self._filters = filters
        self._groups = groups
        self._sorts = sorts
        self._query = query or PaginatedQuery()
        self._mode = mode
        self._page_info: Dict[str, Any] = {
            "page": self._query.page,
            "size": self._query.size,
            "total": 0,
            "pages": 0,
            "list": [],
        }

    def __iter__(self) -> Iterator[Experiment]:
        if self._mode == "get":
            yield from self._iter_paginated()
        else:
            yield from self._iter_filtered()

    def _iter_filtered(self) -> Iterator[Experiment]:
        """POST /runs/shows 模式：复杂过滤，不支持分页。"""
        resp = self._post(
            f"/project/{self._proj_path}/runs/shows",
            data={
                "filters": validate_update_active(self._filters, validate_filter, label="filters"),
                "groups": validate_update_active(self._groups, validate_group, label="groups"),
                "sorts": validate_update_active(self._sorts, validate_sort, label="sorts"),
            },
        )
        if not resp.ok:
            return
        body = resp.data
        runs: list = []
        if isinstance(body, list):
            runs = body
        elif isinstance(body, dict):
            runs = _flatten_runs(body)

        total = len(runs)
        self._page_info.update({"total": total, "page": 1, "size": total})

        for run_data in runs:
            cuid = run_data.get("cuid", "")
            full_path = f"{self._proj_path}/{cuid}"
            yield Experiment(self._ctx, path=full_path, data=run_data)

    def _iter_paginated(self) -> Iterator[Experiment]:
        """GET /runs 模式：标准分页，返回精简信息。"""
        for item in self._paginate(
            f"/project/{self._proj_path}/runs",
            self._query,
            page_info=self._page_info,
        ):
            cuid = item.get("cuid", "")
            full_path = f"{self._proj_path}/{cuid}"
            yield Experiment(
                self._ctx,
                path=full_path,
                data=cast(ApiExperimentType, item),
            )

    def json(self) -> Dict[str, Any]:
        self._page_info["list"] = [r.json() for r in self]
        return self._page_info
