"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Union, cast

from typing_extensions import deprecated

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiResponseType
from swanlab.api.typings.common import (
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricKeyClassLiteral,
    ApiMetricKeyTypeLiteral,
    ApiMetricLogLevelLiteral,
    PaginatedQuery,
    RangeQuery,
)
from swanlab.api.typings.experiment import (
    ApiExperimentLabelType,
    ApiExperimentProfileType,
    ApiExperimentType,
)
from swanlab.api.typings.user import ApiUserType
from swanlab.api.utils import (
    get_properties,
    parse_timestamp_ms,
    resolve_run_path,
    validate_filter,
    validate_group,
    validate_sort,
    validate_update_active,
)
from swanlab.sdk.internal.pkg import console

if TYPE_CHECKING:
    from swanlab.api.series import Series


class Experiment(BaseEntity):
    """
    Represents a SwanLab experiment with lazy-loaded attributes and query / mutation methods.

    Construct with ``data`` to eagerly cache experiment metadata, or ``data=None`` to lazily
    load on first property access. Obtain instances via ``swanlab.Api``::

        exp = api.run("username/project/run_id")
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
        self._cuid: str = (data or {}).get("cuid", "")
        self._data: Optional[ApiExperimentType] = data
        self._project_id = ""

    def _refresh_cuid(self) -> None:
        if self._data:
            self._cuid = self._data.get("cuid", "") or self._cuid

    def _ensure_data(self) -> ApiExperimentType:
        if self._data is None:
            resp = self._get(f"/project/{self._proj_path}/runs/{self._run_slug}")
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
    def root_pro_id(self) -> str:
        return self._ensure_data().get("rootProId", "")

    @property
    def root_exp_id(self) -> str:
        return self._ensure_data().get("rootExpId", "")

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

    @deprecated("Use `series()` method instead.")
    def column(
        self,
        key: str,
        column_class: Optional[ApiColumnClassLiteral] = "CUSTOM",
        column_type: Optional[ApiColumnDataTypeLiteral] = "FLOAT",
    ):
        """
        Get a single column by key under this experiment (fuzzy search, first match).

        Performs fuzzy search (``contains`` on ``name``) and returns the first matching
        column. If multiple columns share a similar name, the first one (ordered by
        ``id DESC``) is returned.

        :param key: Column key, e.g. ``"loss"``, ``"acc"``
        :param column_class: Column class, ``CUSTOM`` (default) or ``SYSTEM``
        :param column_type: Column data type, e.g. ``FLOAT``, ``STRING``, ``IMAGE``
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
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        )

    def metrics(
        self,
        keys: List[str],
        sample: int = 1500,
        ignore_timestamp: bool = False,
        all: bool = False,
        range_query: Optional[Union[Dict[str, Any], RangeQuery]] = None,
    ) -> Dict[str, Any]:
        """
        Fetch scalar metrics (e.g. loss, acc) with three query modes:

        1. **Sampled** (default) — server-side LTTB downsampling, up to ``sample`` data points
        2. **Full** — ``all=True``, no sampling limit
        3. **Range** — filter by step / timestamp / recent time window via ``range_query``

        .. note::
           Modes 2 and 3 (``all`` / ``range_query``) download full-resolution CSV data and
           perform range filtering client-side. Each metric point contains ``step``, ``value``,
           and ``timestamp`` (if available).

        :param keys: Metric keys to fetch, e.g. ``["loss", "acc"]``
        :param sample: Max sampled data points (default 1500, max 1500). Ignored when ``all`` or ``range_query`` is set.
        :param ignore_timestamp: If True, omit timestamp fields from the response
        :param all: If True, fetch full-resolution data without sampling limit
        :param range_query: Range filter — accepts a ``RangeQuery`` object or a plain dict.
            Only supported for SCALAR metrics.

        ---

        **range_query fields**

        ========== ============ =================================================
        Field      Type         Description
        ========== ============ =================================================
        ``type``   ``str``      Filter axis: ``"step"`` (default) or ``"timestamp"``
        ``start``  ``int``      Lower bound (inclusive); None = from beginning
        ``end``    ``int``      Upper bound (inclusive); None = to end
        ``last``   ``int``      Last N milliseconds (mutually exclusive with start/end)
        ``head``   ``int``      First N data points (mutually exclusive with tail)
        ``tail``   ``int``      Last N data points (mutually exclusive with head)
        ========== ============ =================================================

        **Mutual exclusivity**

        - ``head`` and ``tail`` are mutually exclusive
        - ``last`` is mutually exclusive with ``start`` / ``end``
        - ``head`` / ``tail`` can be combined with ``last`` or ``start`` / ``end``

        ---

        **Examples — progressive**

        1. Default sampled query::

               exp.metrics(keys=["loss", "acc"])

        2. Filter by step range::

               exp.metrics(keys=["loss"], range_query={"start": 100, "end": 500})

        3. Step range + first 50 points::

               exp.metrics(keys=["loss"], range_query={"start": 0, "end": 500, "head": 50})

        4. Filter by timestamp (Unix ms, auto-padded if < 13 digits)::

               exp.metrics(keys=["loss"], range_query={
                   "type": "timestamp",
                   "start": 1715769600000,
                   "end": 1715773200000,
               })

        5. Last 5 minutes::

               exp.metrics(keys=["loss"], range_query={"last": 300_000})

        6. Last 5 minutes + first 20 points::

               exp.metrics(keys=["loss"], range_query={"last": 300_000, "head": 20})

        7. Last 30 data points::

               exp.metrics(keys=["loss"], range_query={"tail": 30})
        """
        run_id = self.run_id
        project_id = self.project_id
        from swanlab.api.metric import Metrics

        if isinstance(range_query, dict):
            if range_query.get("type") == "timestamp":
                for key in ("start", "end"):
                    if key in range_query and range_query[key] is not None:
                        range_query[key] = parse_timestamp_ms(range_query[key])
            rq = RangeQuery(**range_query)
        else:
            rq = range_query

        return Metrics(
            ctx=self._ctx,
            project_id=project_id,
            run_id=run_id,
            keys=keys,
            sample=sample,
            metric_type="SCALAR",
            ignore_timestamp=ignore_timestamp,
            all=all,
            range_query=rq,
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
            experiment_name=self.name,
        ).json()

    def summary(
        self,
        keys: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Get summary statistics for scalar metrics (latest, min, max, avg, median, etc.).

        :param keys: Scalar keys to query; None means all keys

        ::

            exp.summary()                 # all keys
            exp.summary(keys=["loss"])    # specific keys
        """
        run_id = self.run_id
        project_id = self.project_id
        from swanlab.api.summary import Summary

        return Summary(
            ctx=self._ctx,
            project_id=project_id,
            experiment_id=run_id,
            keys=keys,
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        ).json()

    def medias(
        self,
        keys: List[str],
        step: Optional[int] = 0,
        all: bool = False,
    ) -> Dict[str, Any]:
        run_id = self.run_id
        project_id = self.project_id
        from swanlab.api.metric import Metrics

        return Metrics(
            ctx=self._ctx,
            project_id=project_id,
            run_id=run_id,
            keys=keys,
            metric_type="MEDIA",
            media_step=step,
            all=all,
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        ).json()

    def logs(
        self,
        offset: Optional[int] = 0,
        level: ApiMetricLogLevelLiteral = "INFO",
        ignore_timestamp: bool = False,
    ) -> Dict[str, Any]:
        run_id = self.run_id
        project_id = self.project_id
        from swanlab.api.metric import Metric

        logs = Metric(
            ctx=self._ctx,
            project_id=project_id,
            run_id=run_id,
            key="LOG",
            log_offset=offset,
            log_level=level,
            metric_type="LOG",
            ignore_timestamp=ignore_timestamp,
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        )
        return logs.json()

    def export_logs(
        self,
        start: int = 0,
        rows: int = 500_000,
    ) -> ApiResponseType:
        """
        Export experiment logs as a .log file.

        :param start: Export start row (0-based), default 0
        :param rows: Number of rows to export, default 500_000, max 500_000
        """
        run_id = self.run_id
        project_id = self.project_id
        from swanlab.api.metric import Metric

        metric = Metric(
            ctx=self._ctx,
            project_id=project_id,
            run_id=run_id,
            key="LOG",
            metric_type="LOG",
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        )
        return metric.export_logs(start=start, rows=rows)

    @deprecated("Use `series()` method instead.")
    def columns(
        self,
        page: int = 1,
        size: int = 100,
        search: Optional[str] = None,
        column_type: Optional[ApiColumnDataTypeLiteral] = None,
        column_class: Optional[ApiColumnClassLiteral] = None,
        all: bool = False,
    ):
        """
        List columns under this experiment (paginated, with optional fuzzy search).

        :param page: Page number, default 1
        :param size: Page size, default 100
        :param search: Fuzzy search keyword (matches column **name**, not key)
        :param column_type: Column type filter, e.g. FLOAT, STRING, IMAGE
        :param column_class: Column class filter, CUSTOM or SYSTEM
        :param all: If True, fetch all pages, default False
        """
        self._ensure_data()
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
            root_pro_id=self.root_pro_id,
            root_exp_id=self.root_exp_id,
        )

    def series(
        self,
        metric_type: ApiMetricKeyTypeLiteral = "SCALAR",
        metric_class: ApiMetricKeyClassLiteral = "CUSTOM",
        search: str = "",
    ) -> "Series":
        """
        List all metric keys for this experiment.

        :param metric_type: ``"SCALAR"`` (default) or ``"MEDIA"``.
        :param metric_class: ``"CUSTOM"`` (default) or ``"SYSTEM"`` — filter keys by class.
        :param search: Fuzzy search filter — case-insensitive substring match on key names.
        :returns: :class:`Series`
        """
        from swanlab.api.series import Series

        run_id = self.run_id
        project_id = self.project_id
        root_pro_id = self.root_pro_id
        root_exp_id = self.root_exp_id

        # 克隆实验的数据存在根实验下，因此查询时直接用根实验的 ID。
        query_pro_id = root_pro_id or project_id
        query_exp_id = root_exp_id or run_id
        experiments: List[Dict[str, str]] = [{"projectId": query_pro_id, "experimentId": query_exp_id}]
        return Series(
            self._ctx,
            experiments=experiments,
            metric_type=metric_type,
            metric_class=metric_class,
            search=search,
            project_id=project_id,
            run_id=run_id,
            root_pro_id=root_pro_id,
            root_exp_id=root_exp_id,
            experiment_name_getter=lambda: self.name,
        )

    def delete(self, commit: bool = False) -> bool:
        """Delete this experiment. ``commit=False`` prints the pending deletion; ``commit=True`` executes it."""
        if not commit:
            name = self.name
            if self._errors:
                return False
            console.warning(f"Experiment to be deleted: run_id: {self._run_slug}, name: {name}")
            return True
        resp = self._delete(f"/project/{self._proj_path}/runs/{self.run_id}")
        return resp.ok

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


def _flatten_runs(runs: Union[list, Dict]) -> list:
    """Flatten grouped experiment data into a flat list of experiments."""
    if isinstance(runs, dict):
        return [item for v in runs.values() for item in _flatten_runs(v)]
    if isinstance(runs, list):
        return list(runs)
    return [runs]


class Experiments(BaseEntity):
    """
    Iterator over experiments in a project.

    Supports two modes:
    - POST mode (default): fetches via ``/runs/shows``, supports complex filtering, no pagination
    - GET mode: fetches via ``/runs``, supports standard pagination, returns summary info

    Usage::

        # POST with complex filters
        for run in api.runs(path="username/project"):
            print(run.name)

        # GET with pagination
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
        """POST /runs/shows mode: complex filtering, no pagination."""
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
        """GET /runs mode: standard pagination, returns summary info."""
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
