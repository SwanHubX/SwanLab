"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple, Union, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import PaginatedQuery
from swanlab.api.typings.experiment import ApiExperimentLabelType, ApiExperimentProfileType, ApiExperimentType
from swanlab.api.typings.user import ApiUserType
from swanlab.api.utils import get_properties, parse_filter


def _resovle_path(path: str) -> Tuple[str, str]:
    """ "path like: user/proj_name/run_id"""
    proj_path, cuid = "", ""
    parts = path.split("/")
    if len(parts) != 3:
        return proj_path, cuid
    cuid = parts[-1]
    proj_path = path.rsplit("/", 1)[0]
    return (
        proj_path,
        cuid,
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
        self._proj_path, self._cuid = _resovle_path(path=path)
        self._data = data

    def _ensure_data(self) -> ApiExperimentType:
        if self._data is None:
            resp = self._get(f"/project/{self._proj_path}/runs/{self._cuid}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiExperimentType, {})
            if not self._cuid and self._data:
                self._cuid = self._data.get("cuid", "")
        return self._data

    @property
    def run_id(self) -> str:
        if self._cuid:
            return self._cuid
        return self._ensure_data().get("cuid", "")

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
        return self._build_web_url(f"@{self._proj_path}/runs/{self.run_id}/chart")

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
        return user_data if isinstance(user_data, dict) else cast(ApiUserType, {})

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
                data = self._data
        return ApiExperimentProfileType(self._ensure_data().get("profile", {}))

    def metrics(
        self, keys: Optional[List[str]] = None, x_axis: Optional[str] = None, sample: Optional[int] = None
    ) -> Any:
        """
        获取实验指标数据，返回 pandas DataFrame。

        :param keys: 指标 key 列表
        :param x_axis: x 轴指标，默认 step
        :param sample: 均匀采样 N 条数据（等间距采样，保留整体趋势）
        """
        from swanlab.vendor import pd

        if not keys:
            return pd.DataFrame()

        fetch_keys = list(keys)
        use_x_axis = x_axis is not None and x_axis != "step"
        if use_x_axis and x_axis is not None:
            fetch_keys.append(x_axis)

        dfs = []
        prefix = ""
        for idx, key in enumerate(fetch_keys):
            resp = self._get(f"/experiment/{self.run_id}/column/csv", params={"key": key})
            if not resp.ok:
                continue
            data = resp.data
            csv_url = data[0].get("url", "") if isinstance(data, list) and data else ""
            if not csv_url:
                continue
            df = pd.read_csv(csv_url, index_col=0)

            if idx == 0:
                first_col = str(df.columns[0])
                suffix = f"{key}_"
                prefix = first_col.split(suffix)[0] if suffix in first_col else ""

            def strip_suffix(col, suffix="_step"):
                return col[: -len(suffix)] if col.endswith(suffix) else col

            df.columns = [
                strip_suffix(col[len(prefix) :]) if prefix and col.startswith(prefix) else strip_suffix(col)
                for col in df.columns
            ]
            dfs.append(df)

        if not dfs:
            return pd.DataFrame()

        result_df = dfs[0].join(dfs[1:], how="outer") if len(dfs) > 1 else dfs[0]
        result_df = result_df.sort_index()

        if use_x_axis:
            result_df = result_df.drop(
                columns=[c for c in result_df.columns if c.endswith("_timestamp")], errors="ignore"
            )
            if x_axis not in result_df.columns:
                return pd.DataFrame()
            cols = [x_axis] + [c for c in result_df.columns if c != x_axis]
            result_df = result_df[cols].dropna(subset=[x_axis])

        if sample is not None and len(result_df) > sample:
            indices = [int(i * (len(result_df) - 1) / (sample - 1)) for i in range(sample)]
            result_df = result_df.iloc[indices]

        return result_df

    def delete(self) -> bool:
        """删除此实验。"""
        resp = self._delete(f"/project/{self._proj_path}/runs/{self._cuid}")
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
        filters: Optional[Dict[str, object]] = None,
        query: Optional[PaginatedQuery] = None,
        mode: str = "post",
    ) -> None:
        super().__init__(ctx)
        self._proj_path = path
        self._filters = filters
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
        parsed_filters = [parse_filter(k, v) for k, v in self._filters.items()] if self._filters else []
        resp = self._post(
            f"/project/{self._proj_path}/runs/shows", data={"filters": parsed_filters, "groups": [], "shows": []}
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
