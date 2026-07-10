"""
@author: caddiesnew
@file: series.py
@time: 2026/7/9
@description: Series 实体类 — 实验指标序列的查询与操作（House key 模型）
"""

from typing import Any, Callable, Dict, Iterator, List, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiMetricKeyClassLiteral, ApiMetricKeyTypeLiteral, ApiResponseType
from swanlab.api.typings.key import ApiSeriesKeyItem
from swanlab.api.utils import get_properties

# 系统指标 key 前缀：SCALAR 类型且以此前缀开头的 key 分类为 SYSTEM
_SYSTEM_KEY_PREFIX = "__swanlab__"


class Key(BaseEntity):
    """A single metric key in a SwanLab experiment — a lightweight handle for querying metric series.

    Unlike Column, Key does not carry rich metadata (name/class/type/createdAt).
    ``metric_type`` is set at construction time and maps directly to the House route
    (scalar / media).

    Usage::

        key = experiment.key("loss", metric_type="SCALAR")
        print(key.json())          # basic info
        data = key.metric()        # scalar line data
        url = key.export_csv()     # CSV export via House
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        project_id: str,
        run_id: str,
        key: str,
        metric_type: ApiMetricKeyTypeLiteral = "SCALAR",
        experiment_name_getter: Optional[Callable[[], str]] = None,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> None:
        super().__init__(ctx)
        self._project_id = project_id
        self._run_id = run_id
        self._key = key
        self._metric_type = metric_type
        self._experiment_name_getter = experiment_name_getter
        self._cached_experiment_name: Optional[str] = None
        self._root_pro_id = root_pro_id
        self._root_exp_id = root_exp_id

    def _resolve_experiment_name(self) -> str:
        if self._cached_experiment_name is not None:
            return self._cached_experiment_name
        if self._experiment_name_getter is not None:
            self._cached_experiment_name = self._experiment_name_getter()
        if self._cached_experiment_name is None:
            self._cached_experiment_name = ""
        return self._cached_experiment_name

    @property
    def project_id(self) -> str:
        return self._project_id

    @property
    def run_id(self) -> str:
        return self._run_id

    @property
    def key(self) -> str:
        return self._key

    @property
    def metric_type(self) -> str:
        return self._metric_type

    @property
    def key_class(self) -> ApiMetricKeyClassLiteral:
        """CUSTOM or SYSTEM — SYSTEM if SCALAR and key starts with __swanlab__."""
        if self._metric_type == "SCALAR" and self._key.startswith(_SYSTEM_KEY_PREFIX):
            return "SYSTEM"
        return "CUSTOM"

    def metric(
        self,
        sample: int = 1500,
        ignore_timestamp: bool = False,
        media_step: Optional[int] = None,
        all: bool = False,
    ) -> Dict[str, Any]:
        from swanlab.api.metric import Metric

        metric = Metric(
            ctx=self._ctx,
            project_id=self._project_id,
            run_id=self._run_id,
            key=self._key,
            sample=sample,
            metric_type=self._metric_type,
            ignore_timestamp=ignore_timestamp,
            media_step=media_step,
            all=all,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        return metric.json()

    def export_csv(self) -> ApiResponseType:
        """Export key data as CSV via House.

        POST /house/metrics/scalar/export → cosKey → presigned download URL.
        Only supported for SCALAR metric_type.
        experiment_name is lazily resolved on first call.
        """
        from swanlab.api.metric import Metric

        if self._metric_type != "SCALAR":
            return ApiResponseType(ok=False, errmsg="export_csv() only support SCALAR metric_type", data=None)

        experiment_name = self._resolve_experiment_name()

        export_column: Dict[str, str] = {
            "experimentId": self._run_id,
            "key": self._key,
        }
        if experiment_name:
            export_column["experimentName"] = experiment_name
        if self._root_pro_id:
            export_column["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            export_column["rootExpId"] = self._root_exp_id

        resp = self._post(
            "/house/metrics/scalar/export",
            data={"projectId": self._project_id, "columns": [export_column]},
        )
        if not resp.ok or not resp.data:
            return resp

        cos_key = resp.data.get("cosKey", "")
        if not cos_key:
            return ApiResponseType(ok=False, errmsg="Invalid response format: missing cosKey", data=None)

        url_map = Metric._fetch_file_presigned_urls(self, [cos_key])
        url = url_map.get(cos_key, "")
        if not url:
            return ApiResponseType(ok=False, errmsg="Failed to get presigned download URL", data=None)
        return ApiResponseType(ok=True, data={"url": url})

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Series(BaseEntity):
    """Metric key listing with full-fetch-and-cache strategy.

    Fetches ALL keys from House on first access, caches them as a flat ``List[str]``,
    then post-filters by ``metric_class``. No pagination state exposed.

    Usage::

        for key in experiment.series():
            print(key.json())

        availability = experiment.series().availability(["loss", "acc"])
    """

    _PAGE_SIZE = 2000

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        experiments: List[Dict[str, str]],
        metric_type: ApiMetricKeyTypeLiteral = "SCALAR",
        metric_class: ApiMetricKeyClassLiteral = "CUSTOM",
        search: str = "",
        project_id: str = "",
        run_id: str = "",
        root_pro_id: str = "",
        root_exp_id: str = "",
        experiment_name_getter: Optional[Callable[[], str]] = None,
    ) -> None:
        super().__init__(ctx)
        if metric_type not in ("SCALAR", "MEDIA"):
            raise ValueError(f"Invalid metric_type: {metric_type!r}, expected 'SCALAR' or 'MEDIA'")
        if metric_class not in ("CUSTOM", "SYSTEM"):
            raise ValueError(f"Invalid metric_class: {metric_class!r}, expected 'CUSTOM' or 'SYSTEM'")
        if not isinstance(experiments, list) or not experiments:
            raise ValueError("experiments must be a non-empty list of dicts")

        self._experiments = experiments
        self._metric_type: ApiMetricKeyTypeLiteral = metric_type
        self._metric_class: ApiMetricKeyClassLiteral = metric_class
        self._search = search
        self._project_id = project_id
        self._run_id = run_id
        self._root_pro_id = root_pro_id
        self._root_exp_id = root_exp_id
        self._experiment_name_getter = experiment_name_getter
        self._all_keys: Optional[List[str]] = None
        self._cached_list: Optional[List[Key]] = None

    # ------------------------------------------------------------------
    # 内部：全量拉取 + 缓存
    # ------------------------------------------------------------------

    def _fetch_all_keys(self) -> List[str]:
        """Fetch all keys across all pages (cursor pagination), cache and return."""
        if self._all_keys is not None:
            return self._all_keys
        keys: List[str] = []
        cursor = ""
        path = f"/house/metrics/{self._metric_type.lower()}/keys"
        while True:
            resp = self._post(
                path,
                data={
                    "experiments": self._experiments,
                    "limit": self._PAGE_SIZE,
                    "cursor": cursor,
                },
            )
            if not resp.ok or not resp.data:
                break
            page: Dict[str, Any] = resp.data
            keys.extend(page.get("keys", []))
            if not page.get("hasMore", False) or not page.get("nextCursor", ""):
                break
            cursor = page["nextCursor"]
        self._all_keys = keys
        return keys

    def _filtered_keys(self) -> List[str]:
        """Apply metric_class + search post-filter on cached full key list."""
        keys = self._fetch_all_keys()
        want_system = self._metric_class == "SYSTEM"
        if self._metric_type == "SCALAR":
            result = [k for k in keys if k.startswith(_SYSTEM_KEY_PREFIX) == want_system]
        else:
            # MEDIA: 全为 CUSTOM
            result = [] if want_system else keys
        # 模糊搜索：大小写不敏感的子串匹配
        if self._search:
            needle = self._search.lower()
            result = [k for k in result if needle in k.lower()]
        return result

    def _ensure_batch(self) -> List["Key"]:
        if self._cached_list is not None:
            return self._cached_list
        self._cached_list = [
            Key(
                self._ctx,
                project_id=self._project_id,
                run_id=self._run_id,
                key=key_str,
                metric_type=self._metric_type,
                experiment_name_getter=self._experiment_name_getter,
                root_pro_id=self._root_pro_id,
                root_exp_id=self._root_exp_id,
            )
            for key_str in self._filtered_keys()
        ]
        return self._cached_list

    # ------------------------------------------------------------------
    # 公开接口
    # ------------------------------------------------------------------

    def __iter__(self) -> Iterator[Key]:
        yield from self._ensure_batch()

    def availability(self, keys: List[str]) -> Dict[str, List[str]]:
        """Check which experiments have data for the given keys.

        POST /house/metrics/{scalar|media}/availability → {key: experimentId[]}.
        An empty array means no experiment has this key.
        """
        if not isinstance(keys, list) or not keys:
            return {}
        path = f"/house/metrics/{self._metric_type.lower()}/availability"
        resp = self._post(path, data={"keys": keys, "experiments": self._experiments})
        if resp.ok and isinstance(resp.data, dict):
            return resp.data
        return {}

    @property
    def total(self) -> int:
        """Total number of keys after metric_class filtering (triggers fetch)."""
        return len(self._filtered_keys())

    def json(self) -> Dict[str, Any]:
        keys = self._filtered_keys()
        if self._metric_type == "SCALAR":
            result_list: List[ApiSeriesKeyItem] = [
                {"key": k, "metric_class": "SYSTEM" if k.startswith(_SYSTEM_KEY_PREFIX) else "CUSTOM"} for k in keys
            ]
        else:
            result_list = [{"key": k, "metric_class": "CUSTOM"} for k in keys]
        return {
            "list": result_list,
            "metricType": self._metric_type,
            "projectId": self._project_id,
            "runId": self._run_id,
        }
