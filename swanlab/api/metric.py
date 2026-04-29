"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/20
@description: Metric 实体类 — 指标序列的查询与操作
"""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiColumnCsvExportType, ApiResponseType
from swanlab.api.typings.common import ApiMetricColumnTypeLiteral, ApiMetricLogLevelLiteral
from swanlab.api.typings.metric import (
    ApiLogSeriesType,
    ApiMediaItemDataType,
    ApiMediaSeriesType,
    ApiScalarSeriesType,
)
from swanlab.api.utils import get_properties, validate_metric_keys, validate_metric_log_level, validate_metric_type
from swanlab.sdk.internal.pkg import console

_SCALAR_STATISTIC_FIELDS = ("min", "max", "avg", "median", "latest")
_METRIC_SHARED_KEYS = frozenset({"project_id", "run_id", "metric_type"})


def _extract_csv_url(data: Any) -> str:
    if isinstance(data, list) and data:
        return data[0].get("url", "")
    if isinstance(data, dict):
        return data.get("url", "")
    return ""


class Metric(BaseEntity):
    """
    表示一个 SwanLab 指标列（非单个数值，而是一组序列）。

    支持 SCALAR / MEDIA / LOG 三种类型，按需 Lazy Loading。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        project_id: str,
        run_id: str,
        key: Optional[str] = "",
        sample: int = 1000,
        log_offset: Optional[int] = 0,  # 标记第几个分片，仅对 Log metric_type 有效
        log_level: ApiMetricLogLevelLiteral = "INFO",
        metric_type: str = "SCALAR",
        data: Optional[Dict[str, Any]] = None,
        ignore_timestamp: bool = False,
        media_step: Optional[int] = None,
        all: bool = False,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> None:
        super().__init__(ctx)
        validate_metric_type(metric_type, key)
        if metric_type == "LOG":
            validate_metric_log_level(log_level)
        self._project_id = project_id
        self._run_id = run_id
        self._key = key
        self._data: Optional[Dict[str, Any]] = data
        self._metric_type = metric_type
        self._ignore_timestamp = ignore_timestamp
        # 采样值， scalar 时生效
        self._sample = sample
        # 偏移量，仅对 Log metric_type 有效， 默认为 0
        self._offset = log_offset
        self._log_level = log_level
        self._media_step = media_step
        self._all = all
        self._root_pro_id = root_pro_id
        self._root_exp_id = root_exp_id

    # 类型 → 加载方法 的分发表，新增类型只需在此注册
    _FETCH_DISPATCH = {
        "SCALAR": "_fetch_scalar",
        "MEDIA": "_fetch_media",
        "LOG": "_fetch_logs",
    }

    def _ensure_data(self) -> Dict[str, Any]:
        if self._data is None:
            if self._metric_type == "MEDIA" and self._all:
                method_name = "_fetch_media_all"
            else:
                method_name = self._FETCH_DISPATCH.get(self._metric_type, "_fetch_scalar")
            self._data = getattr(self, method_name)()
        assert self._data is not None
        return self._data

    @property
    def project_id(self) -> str:
        return self._project_id

    @property
    def run_id(self) -> str:
        return self._run_id

    @property
    def key(self) -> str:
        return self._key or ""

    @property
    def metric_type(self) -> str:
        return self._metric_type

    @property
    def metrics(self) -> List[Any]:
        return self._ensure_data().get("metrics", [])

    @property
    def logs(self) -> List[Any]:
        return self._ensure_data().get("logs", [])

    @property
    def count(self) -> int:
        return self._ensure_data().get("count", 0)

    @property
    def steps(self) -> List[int]:
        return self._ensure_data().get("steps", [])

    # only available for media type
    @property
    def step(self) -> Optional[int]:
        return self._ensure_data().get("step")

    # ------------------------------------------------------------------
    # 请求辅助函数
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_first(resp: ApiResponseType) -> Optional[Dict[str, Any]]:
        """从列表型 API 响应中提取第一个元素，失败返回 None。"""
        if resp.ok and isinstance(resp.data, list) and resp.data:
            return resp.data[0]
        return None

    @staticmethod
    def _build_column_ref(
        experiment_id: str,
        key: str,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> Dict[str, str]:
        ref: Dict[str, str] = {"experimentId": experiment_id, "key": key}
        if root_pro_id:
            ref["rootProId"] = root_pro_id
        if root_exp_id:
            ref["rootExpId"] = root_exp_id
        return ref

    @staticmethod
    def _build_scalar_payload(
        project_id: str,
        run_id: str,
        keys: List[str],
        sample: int = 1500,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> Dict[str, Any]:
        return {
            "projectId": project_id,
            "xType": "step",
            "range": [0, 0],
            "columns": [Metric._build_column_ref(run_id, key, root_pro_id, root_exp_id) for key in keys],
            "num": sample if sample <= 1500 else 1500,
        }

    @staticmethod
    def _build_media_payload(
        project_id: str,
        run_id: str,
        keys: List[str],
        step: Optional[int] = None,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "projectId": project_id,
            "columns": [Metric._build_column_ref(run_id, key, root_pro_id, root_exp_id) for key in keys],
        }
        if step is not None:
            payload["step"] = step
        return payload

    def _build_log_params(self) -> Dict[str, Any]:
        params: Dict[str, Any] = {
            "projectId": self.project_id,
            "experimentId": self.run_id,
            "size": 1000,
            "epoch": self._offset,
            "level": self._log_level,
        }
        if self._root_pro_id:
            params["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            params["rootExpId"] = self._root_exp_id
        return params

    # ------------------------------------------------------------------
    # 类型专属加载
    # ------------------------------------------------------------------

    def _fetch_scalar(self) -> ApiScalarSeriesType:
        res = ApiScalarSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        if self._root_pro_id:
            res["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            res["rootExpId"] = self._root_exp_id
        payload = self._build_scalar_payload(
            self.project_id,
            self.run_id,
            [self.key],
            self._sample,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )

        # 1. 获取折线数据
        raw_data = self._extract_first(self._post("/house/metrics/scalar", data=payload))
        if raw_data is None:
            return res
        res["metrics"] = raw_data.get("metrics", [])

        # 2. 获取统计值
        stat_data = self._extract_first(self._post("/house/metrics/scalar/value", data=payload))
        if stat_data is None:
            return res
        for field in _SCALAR_STATISTIC_FIELDS:
            res[field] = stat_data.get(field, {})
        return res

    @staticmethod
    def _fetch_presigned_urls(entity: BaseEntity, prefix: str, paths: List[str]) -> Dict[str, str]:
        """批量获取预签名下载链接，返回 path → url 映射。"""
        if not paths:
            return {}
        resp = entity._post("/resources/presigned/get", data={"prefix": prefix, "paths": paths})
        if not resp.ok or not isinstance(resp.data, dict):
            return {}
        urls = resp.data.get("urls", [])
        return dict(zip(paths, urls)) if urls else {}

    @staticmethod
    def _build_media_items(
        entry: Dict[str, Any],
        url_map: Dict[str, str],
    ) -> List[ApiMediaItemDataType]:
        """将单个 metric entry 的 data/more 合并为 items，注入预签名 url。"""
        paths = entry.get("data", [])
        mores = entry.get("more", [])
        items: List[ApiMediaItemDataType] = []
        for i, path in enumerate(paths):
            item: ApiMediaItemDataType = {}
            if path in url_map:
                item["url"] = url_map[path]
            if i < len(mores) and isinstance(mores[i], dict):
                item.update(mores[i])
            items.append(item)
        return items

    def _fetch_media(self) -> ApiMediaSeriesType:
        res = ApiMediaSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = self._build_media_payload(
            self.project_id,
            self.run_id,
            [self.key],
            step=self._media_step,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/media", data=payload)
        if not raw_resp.ok or not raw_resp.data:
            return res
        data = raw_resp.data
        if not isinstance(data, dict):
            return res

        res["steps"] = data.get("steps", [])
        step_val = data.get("step")
        if step_val is not None:
            res["step"] = step_val

        metrics_raw: List[Dict[str, Any]] = data.get("metrics", [])
        metric_entry = next((m for m in metrics_raw if m.get("key") == self.key), None)
        if metric_entry is None:
            return res

        prefix = f"{self.project_id}/{self.run_id}"
        all_paths = metric_entry.get("data", [])
        url_map = self._fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched: run_id[{self.run_id}], key[{self.key}] - {len(all_paths)} items, requesting presigned urls..."
            )
        items = self._build_media_items(metric_entry, url_map)
        res["metrics"] = [{"index": data.get("step", 0), "items": items}]
        return res

    def _fetch_media_all(self) -> ApiMediaSeriesType:
        res = ApiMediaSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = self._build_media_payload(
            self.project_id,
            self.run_id,
            [self.key],
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/f_media", data=payload)
        raw_data = self._extract_first(raw_resp)
        if raw_data is None:
            return res

        prefix = f"{self.project_id}/{self.run_id}"
        all_paths = [p for entry in raw_data.get("metrics", []) for p in entry.get("data", [])]
        url_map = self._fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched (all): run_id[{self.run_id}], key[{self.key}] - {len(all_paths)} items, requesting presigned urls..."
            )
        res["metrics"] = [
            {"index": entry.get("index", 0), "items": self._build_media_items(entry, url_map)}
            for entry in raw_data.get("metrics", [])
        ]
        return res

    def _fetch_logs(self) -> ApiLogSeriesType:
        res = ApiLogSeriesType(projectId=self.project_id, experimentId=self.run_id, key="LOG")
        params = self._build_log_params()
        raw_resp = self._get("/house/metrics/log", params=params)
        if not raw_resp.ok or not raw_resp.data:
            return res
        data = raw_resp.data
        if isinstance(data, dict):
            res["logs"] = data.get("logs", [])
            res["count"] = data.get("count", 0)
        return res

    # ------------------------------------------------------------------
    # 导出
    # ------------------------------------------------------------------

    def export_csv(self) -> ApiResponseType:
        """导出列数据为 CSV。"""
        if self.metric_type != "SCALAR":
            return ApiResponseType(ok=False, errmsg="export_csv() only support SCALAR metric_type", data=None)
        resp = self._get(f"/experiment/{self._run_id}/column/csv", params={"key": self.key})
        if not resp.ok:
            return resp
        url = _extract_csv_url(resp.data)
        if not url:
            return ApiResponseType(ok=False, errmsg="Invalid response format", data=None)
        return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))

    def json(self) -> Dict[str, Any]:
        result = get_properties(self)
        data = self._ensure_data()

        if self._metric_type == "SCALAR":
            if "url" in data:
                result.pop("metrics", None)
                result["url"] = data["url"]
            for field in _SCALAR_STATISTIC_FIELDS:
                val = data.get(field)
                if val:
                    result[field] = val

        if self._metric_type == "LOG":
            result.pop("metrics", None)
        else:
            result.pop("logs", None)
            result.pop("count", None)

        if self._metric_type != "MEDIA" or "steps" not in data:
            result.pop("steps", None)
            result.pop("step", None)

        if self._ignore_timestamp:
            items = result.get("metrics", []) or result.get("logs", [])
            if isinstance(items, list):
                for item in items:
                    if isinstance(item, dict):
                        item.pop("timestamp", None)

        return result


class Metrics(BaseEntity):
    """
    批量指标数据的迭代器。

    一次 metrics 查询只支持一种 metric_type（SCALAR 或 MEDIA），不支持 LOG。
    通过 payload 的 columns 数组一次性传递多个 key，减少网络请求。

    用法::

        for m in experiment.metrics(keys=["loss", "acc"], metric_type="SCALAR"):
            print(m.key, m.metrics)
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        project_id: str,
        run_id: str,
        keys: List[str],
        metric_type: ApiMetricColumnTypeLiteral,
        sample: int = 1500,
        ignore_timestamp: bool = False,
        media_step: Optional[int] = None,
        all: bool = False,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> None:
        super().__init__(ctx)
        validate_metric_keys(keys)
        validate_metric_type(metric_type, keys[0])
        if metric_type == "LOG":
            raise ValueError("Metrics does not support LOG metric_type, use Experiment.logs() instead")
        self._project_id = project_id
        self._run_id = run_id
        self._keys = keys
        self._metric_type = metric_type

        self._ignore_timestamp = ignore_timestamp
        self._media_step = media_step
        self._all = all
        self._root_pro_id = root_pro_id
        self._root_exp_id = root_exp_id
        self._page_info: Dict[str, Any] = {
            "keys": keys,
            "metricType": metric_type,
            "projectId": project_id,
            "experimentId": run_id,
            "list": [],
        }
        self._sample = sample
        if sample > 1500:
            console.warning(f"Get sample = [{sample}], expected <= 1500, will be constrainted automatically..")
            self._sample = 1500

    def __iter__(self) -> Iterator[Metric]:
        if self._metric_type == "SCALAR":
            if self._all:
                yield from self._fetch_scalars_all()
            else:
                yield from self._fetch_scalars()
        else:
            if self._all:
                yield from self._fetch_medias_all()
            else:
                yield from self._fetch_medias()

    def _build_metric(self, key: str, data: Dict[str, Any]) -> Metric:
        return Metric(
            ctx=self._ctx,
            project_id=self._project_id,
            run_id=self._run_id,
            key=key,
            metric_type=self._metric_type,
            sample=self._sample,
            ignore_timestamp=self._ignore_timestamp,
            media_step=self._media_step,
            data=data,
            all=self._all,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )

    def _fetch_scalars(self) -> Iterator[Metric]:
        payload = Metric._build_scalar_payload(
            self._project_id,
            self._run_id,
            self._keys,
            self._sample,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )

        # 1. 获取折线数据
        scalar_resp = self._post("/house/metrics/scalar", data=payload)
        scalar_list: List[Dict[str, Any]] = (
            scalar_resp.data if scalar_resp.ok and isinstance(scalar_resp.data, list) else []
        )

        # 2. 获取统计值
        value_resp = self._post("/house/metrics/scalar/value", data=payload)
        value_list: List[Dict[str, Any]] = value_resp.ok and isinstance(value_resp.data, list) and value_resp.data or []

        for i, key in enumerate(self._keys):
            data: Dict[str, Any] = {
                "projectId": self._project_id,
                "experimentId": self._run_id,
                "key": key,
                "metrics": [],
            }
            if i < len(scalar_list):
                data["metrics"] = scalar_list[i].get("metrics", [])
            if i < len(value_list):
                for field in _SCALAR_STATISTIC_FIELDS:
                    val = value_list[i].get(field)
                    if val is not None:
                        data[field] = val
            yield self._build_metric(key, data)

    def _fetch_scalars_all(self) -> Iterator[Metric]:
        urls: Dict[str, str] = {}
        for key in self._keys:
            resp = self._get(f"/experiment/{self._run_id}/column/csv", params={"key": key})
            if resp.ok and resp.data:
                urls[key] = _extract_csv_url(resp.data)

        payload = Metric._build_scalar_payload(
            self._project_id,
            self._run_id,
            self._keys,
            self._sample,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        value_resp = self._post("/house/metrics/scalar/value", data=payload)
        value_list: List[Dict[str, Any]] = value_resp.ok and isinstance(value_resp.data, list) and value_resp.data or []

        for i, key in enumerate(self._keys):
            data: Dict[str, Any] = {
                "projectId": self._project_id,
                "experimentId": self._run_id,
                "key": key,
                "url": urls.get(key, ""),
            }
            if i < len(value_list):
                for field in _SCALAR_STATISTIC_FIELDS:
                    val = value_list[i].get(field)
                    if val is not None:
                        data[field] = val
            yield self._build_metric(key, data)

    def _fetch_medias(self) -> Iterator[Metric]:
        payload = Metric._build_media_payload(
            self._project_id,
            self._run_id,
            self._keys,
            step=self._media_step,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/media", data=payload)
        if not raw_resp.ok or not raw_resp.data:
            return
        resp_data = raw_resp.data
        if not isinstance(resp_data, dict):
            return

        steps = resp_data.get("steps", [])
        current_step = resp_data.get("step")
        metrics_raw: List[Dict[str, Any]] = resp_data.get("metrics", [])

        prefix = f"{self._project_id}/{self._run_id}"
        all_paths = [p for entry in metrics_raw for p in entry.get("data", [])]
        url_map = Metric._fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched: run_id[{self._run_id}] - {len(all_paths)} items across {len(self._keys)} keys, requesting presigned urls..."
            )

        key_to_entry: Dict[str, Dict[str, Any]] = {e.get("key", ""): e for e in metrics_raw}
        for key in self._keys:
            data: Dict[str, Any] = {
                "projectId": self._project_id,
                "experimentId": self._run_id,
                "key": key,
                "steps": steps,
                "step": current_step,
                "metrics": [],
            }
            entry = key_to_entry.get(key)
            if entry:
                items = Metric._build_media_items(entry, url_map)
                data["metrics"] = [{"index": current_step or 0, "items": items}]
            yield self._build_metric(key, data)

    def _fetch_medias_all(self) -> Iterator[Metric]:
        payload = Metric._build_media_payload(
            self._project_id,
            self._run_id,
            self._keys,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/f_media", data=payload)
        if not raw_resp.ok or not raw_resp.data:
            return
        raw_list = raw_resp.data
        if not isinstance(raw_list, list):
            return

        prefix = f"{self._project_id}/{self._run_id}"
        all_paths = [p for entry in raw_list for m in entry.get("metrics", []) for p in m.get("data", [])]
        url_map = Metric._fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched (all): run_id[{self._run_id}] - {len(all_paths)} items across {len(self._keys)} keys, requesting presigned urls..."
            )

        key_to_entry: Dict[str, Dict[str, Any]] = {e.get("key", ""): e for e in raw_list}
        for key in self._keys:
            data: Dict[str, Any] = {
                "projectId": self._project_id,
                "experimentId": self._run_id,
                "key": key,
                "metrics": [],
            }
            entry = key_to_entry.get(key)
            if entry:
                metrics_list: List[Dict[str, Any]] = []
                for m in entry.get("metrics", []):
                    items = Metric._build_media_items(m, url_map)
                    metrics_list.append({"index": m.get("index", 0), "items": items})
                data["metrics"] = metrics_list
            yield self._build_metric(key, data)

    def json(self) -> Dict[str, Any]:
        self._page_info["list"] = [{k: v for k, v in m.json().items() if k not in _METRIC_SHARED_KEYS} for m in self]
        return self._page_info
