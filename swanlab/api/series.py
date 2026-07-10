"""
@author: caddiesnew
@file: series.py
@time: 2026/7/9
@description: Series 实体类 — 实验指标序列查询
"""

from typing import Any, Callable, Dict, Iterator, List, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiMetricKeyClassLiteral, ApiMetricKeyTypeLiteral, ApiResponseType
from swanlab.api.utils import get_properties

# 系统指标 key 前缀：SCALAR 类型且以此前缀开头的 key 分类为 SYSTEM
_SYSTEM_KEY_PREFIX = "__swanlab__"


class Key(BaseEntity):
    """A single metric key — lightweight handle for querying metric data and exporting CSV.

    Usage::

        key = experiment.key("loss")
        print(key.json())          # {key, metric_type, key_class, ...}
        data = key.metric()        # scalar line data
        url = key.export_csv()     # CSV download URL
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
        """Key classification: ``"CUSTOM"`` or ``"SYSTEM"``."""
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
        """Fetch metric data points for this key.

        :param sample: Max data points (downsampled server-side). Default 1500, max 1500.
        :param ignore_timestamp: If True, omit ``timestamp`` field from each data point.
        :param media_step: Step filter for MEDIA metrics. If None, returns all steps.
        :param all: If True, fetch full-resolution data without downsampling.
        :returns: ``{"list": [{"step", "value", "timestamp", "key"}], ...}``
        """
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
        """Export this key's data as a CSV file (SCALAR only).

        :returns: ``ApiResponseType(ok=True, data={"url": "<download_url>"})``.
            Returns ``ok=False`` for MEDIA keys.
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
    """Metric key listing for one or more experiments.

    Iterating yields :class:`Key` objects. Use :meth:`json` for a plain-dict summary.

    Usage::

        series = experiment.series()
        for key in series:              # iterate Key objects
            print(key.json())
        series.availability(["loss"])   # check key existence
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
        self._cached_filtered: Optional[List[str]] = None
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
        """Apply metric_class + search post-filter on cached full key list. Result is cached."""
        if self._cached_filtered is not None:
            return self._cached_filtered
        keys = self._fetch_all_keys()
        want_system = self._metric_class == "SYSTEM"
        # 分支提循环外，用 not 替代 == bool 比较
        if self._metric_type == "SCALAR":
            if want_system:
                result = [k for k in keys if k.startswith(_SYSTEM_KEY_PREFIX)]
            else:
                result = [k for k in keys if not k.startswith(_SYSTEM_KEY_PREFIX)]
        else:
            # MEDIA: 全为 CUSTOM
            result = [] if want_system else keys
        # 模糊搜索：大小写不敏感的子串匹配
        if self._search:
            needle = self._search.lower()
            result = [k for k in result if needle in k.lower()]
        self._cached_filtered = result
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
        """Check which of the given keys exist in the experiment(s).

        :param keys: Key names to check.
        :returns: ``{key_name: [experiment_id, ...]}``. An empty list means no experiment has that key.
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
        """Number of keys after filtering."""
        return len(self._filtered_keys())

    def json(self) -> Dict[str, Any]:
        """Serialize to a JSON-compatible dict.

        :returns:

            ========== ============ ===========================================
            Field       Type         Description
            ========== ============ ===========================================
            ``keys``    ``list[str]`` Filtered key names.
            ``metricClass`` ``str``   ``"CUSTOM"`` or ``"SYSTEM"``.
            ``metricType`` ``str``    ``"SCALAR"`` or ``"MEDIA"``.
            ``projectId`` ``str``     Project ID.
            ``runId``   ``str``       Experiment (run) ID.
            ========== ============ ===========================================
        """
        return {
            "keys": self._filtered_keys(),
            "metricClass": self._metric_class,
            "metricType": self._metric_type,
            "projectId": self._project_id,
            "runId": self._run_id,
        }
