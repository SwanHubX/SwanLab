"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Union, cast

from swanlab.utils import parse_column_type, to_camel_case

from .base import BaseEntity
from .typings.experiment import ApiExperimentType, ApiExperimentUserType
from .utils import Label, get_properties

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Profile:
    """Experiment profile containing config, metadata, requirements, and conda info."""

    def __init__(self, data: Dict) -> None:
        self._data = data

    @staticmethod
    def _clean_field(value: Any) -> Any:
        """Recursively clean config field, removing desc/sort and keeping value."""
        if isinstance(value, dict):
            if "value" in value:
                return Profile._clean_field(value["value"])
            else:
                return {k: Profile._clean_field(v) for k, v in value.items()}
        elif isinstance(value, list):
            return [Profile._clean_field(item) for item in value]
        return value

    @property
    def config(self) -> Dict:
        """Experiment configuration (cleaned, without desc/sort fields)."""
        raw_config = self._data.get("config", {})
        return {k: Profile._clean_field(v) for k, v in raw_config.items()} if isinstance(raw_config, dict) else {}

    @property
    def metadata(self) -> Dict:
        return self._data.get("metadata", {})

    @property
    def requirements(self) -> str:
        return self._data.get("requirements", "")

    @property
    def conda(self) -> str:
        return self._data.get("conda", "")


class Experiment(BaseEntity):
    """
    表示一个 SwanLab 实验。

    支持双模式：构造时传入 data，或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        path: str,
        data: Optional[ApiExperimentType] = None,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._path = path  # 'username/project-name'
        self._data = data

    def _ensure_data(self) -> ApiExperimentType:
        if self._data is None:
            resp = self._get(f"/project/{self._path}/runs/{self.id}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiExperimentType, {})
        return self._data

    @property
    def id(self) -> str:
        return self._data.get("cuid", "") if self._data is not None else self._ensure_data().get("cuid", "")

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def description(self) -> str:
        return self._ensure_data().get("description", "")

    @property
    def state(self) -> str:
        return self._ensure_data().get("state", "")

    @property
    def url(self) -> str:
        return self._build_url(f"@{self._path}/runs/{self.id}/chart")

    @property
    def show(self) -> bool:
        return self._ensure_data().get("show", True)

    @property
    def labels(self) -> List[Label]:
        return [Label(label["name"]) for label in self._ensure_data().get("labels", [])]

    @property
    def group(self) -> str:
        return self._ensure_data().get("cluster", "")

    @property
    def job_type(self) -> str:
        return self._ensure_data().get("job", "")

    @property
    def user(self) -> ApiExperimentUserType:
        user_data = self._ensure_data().get("user", {})
        return user_data if isinstance(user_data, dict) else cast(ApiExperimentUserType, {})

    @property
    def created_at(self) -> str:
        return self._ensure_data().get("createdAt", "")

    @property
    def finished_at(self) -> str:
        return self._ensure_data().get("finishedAt", "")

    @property
    def profile(self) -> Profile:
        """Experiment profile containing config, metadata, requirements, and conda."""
        data = self._ensure_data()
        if "profile" not in data and self.id:
            resp = self._get(f"/project/{self._path}/runs/{self.id}")
            if resp.ok and resp.data:
                self._data = resp.data
                data = self._data
        return Profile(data.get("profile", {}))

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
            resp = self._get(f"/experiment/{self.id}/column/csv", params={"key": key})
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
        resp = self._delete(f"/project/{self._path}/runs/{self.id}")
        return resp.ok

    def to_dict(self) -> Dict[str, Any]:
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

    用法::

        for run in api.runs("username/project"):
            print(run.name)
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        path: str,
        filters: Optional[Dict[str, object]] = None,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._path = path
        self._filters = filters

    def __iter__(self) -> Iterator[Experiment]:
        parsed_filters = (
            [
                {
                    "key": to_camel_case(key) if parse_column_type(key) == "STABLE" else key.split(".", 1)[-1],
                    "active": True,
                    "value": [value],
                    "op": "EQ",
                    "type": parse_column_type(key),
                }
                for key, value in self._filters.items()
            ]
            if self._filters
            else []
        )
        resp = self._post(f"/project/{self._path}/runs/shows", data={"filters": parsed_filters})
        if not resp.ok:
            return
        body = resp.data
        runs: list = []
        if isinstance(body, list):
            runs = body
        elif isinstance(body, dict):
            runs = _flatten_runs(body)

        for run_data in runs:
            yield Experiment(self._client, self._web_host, self._api_host, path=self._path, data=run_data)

    def to_dict(self) -> Dict[str, Any]:
        return {"path": self._path}
