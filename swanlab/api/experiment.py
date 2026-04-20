"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional

from .base import BaseEntity
from .typings.experiment import ApiExperimentType
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
            # path 格式为 'username/project'，需要构造完整请求路径
            self._data = self._get(f"/project/{self._path}/runs/{self.id}" if self.id else f"/project/{self._path}")
        return self._data

    @property
    def id(self) -> str:
        if self._data is not None:
            return self._data.get("cuid", "")
        # data 为空时无法获取 id，触发一次 fetch
        # 先尝试不带 id 的路径来获取
        return self._ensure_data().get("cuid", "")

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
    def user(self) -> Any:
        user_data = self._ensure_data().get("user", {})
        return user_data if isinstance(user_data, dict) else {}

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
        if "profile" not in data:
            # 重新 fetch 获取完整数据
            self._data = self._get(f"/project/{self._path}/runs/{self.id}")
            data = self._data
        return Profile(data.get("profile", {}))

    def metrics(
        self, keys: Optional[List[str]] = None, x_axis: Optional[str] = None, sample: Optional[int] = None
    ) -> Any:
        """
        获取实验指标数据，返回 pandas DataFrame。

        :param keys: 指标 key 列表
        :param x_axis: x 轴指标，默认 step
        :param sample: 返回前 N 条数据
        """
        try:
            import pandas as pd
        except ImportError:
            raise TypeError("pandas is required for metrics(). Install with 'pip install pandas'.")

        if not keys:
            return pd.DataFrame()

        fetch_keys = list(keys)
        use_x_axis = x_axis is not None and x_axis != "step"
        if use_x_axis and x_axis is not None:
            fetch_keys.append(x_axis)

        dfs = []
        prefix = ""
        for idx, key in enumerate(fetch_keys):
            resp = self._client.get(f"/experiment/{self.id}/column/csv", params={"key": key})
            csv_url = resp.data[0].get("url", "") if isinstance(resp.data, list) and resp.data else ""
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
                raise ValueError(f"x_axis '{x_axis}' not found in result DataFrame")
            cols = [x_axis] + [c for c in result_df.columns if c != x_axis]
            result_df = result_df[cols].dropna(subset=[x_axis])

        if sample is not None:
            result_df = result_df.head(sample)

        return result_df

    def delete(self) -> None:
        """删除此实验。"""
        self._delete(f"/project/{self._path}/runs/{self.id}")

    def to_dict(self) -> Dict[str, Any]:
        return get_properties(self)
