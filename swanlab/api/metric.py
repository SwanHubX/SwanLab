"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/20
@description: Metric 实体类 — 指标序列的查询与操作
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, Iterator, List, Optional, Union

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiColumnCsvExportType, ApiResponseType
from swanlab.api.typings.common import (
    MAX_METRIC_KEY_BATCH_SIZE,
    ApiMetricColumnTypeLiteral,
    ApiMetricLogLevelLiteral,
    RangeQuery,
)
from swanlab.api.typings.metric import (
    ApiLogSeriesType,
    ApiMediaItemDataType,
    ApiMediaSeriesType,
    ApiScalarSeriesType,
)
from swanlab.api.utils import get_properties, validate_metric_keys, validate_metric_log_level, validate_metric_type
from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client

_SCALAR_STATISTIC_FIELDS = ("min", "max", "avg", "median", "latest")
_METRIC_SHARED_KEYS = frozenset({"project_id", "run_id", "metric_type"})


def _align_entries_by_key(entries: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """将后端返回的列表按 ``key`` 字段映射为 dict，应对后端可能省略列或乱序返回。"""
    indexed: Dict[str, Dict[str, Any]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        key = entry.get("key")
        if isinstance(key, str) and key:
            indexed[key] = entry
    return indexed


def _merge_value_stats(
    step_list: List[Dict[str, Any]],
    time_list: List[Dict[str, Any]],
    keys: List[str],
) -> List[Dict[str, Any]]:
    """合并 step/time 两种 x_type 的 value stat 响应为 per-key 统计字典。"""
    step_by_key = _align_entries_by_key(step_list)
    time_by_key = _align_entries_by_key(time_list)
    merged: List[Dict[str, Any]] = []
    for key in keys:
        step_entry = step_by_key.get(key, {})
        time_entry = time_by_key.get(key, {})
        entry: Dict[str, Any] = {}
        for field in _SCALAR_STATISTIC_FIELDS:
            step_val = step_entry.get(field)
            time_val = time_entry.get(field)
            if isinstance(step_val, dict):
                stat = dict(step_val)
                if isinstance(time_val, dict) and time_val.get("index") is not None:
                    stat["timestamp"] = time_val["index"]
                entry[field] = stat
            elif isinstance(time_val, dict):
                entry[field] = dict(time_val)
        merged.append(entry)
    return merged


def _extract_csv_url(data: Any) -> str:
    if isinstance(data, list) and data:
        return _extract_csv_url(data[0])
    if isinstance(data, dict):
        url = data.get("url", "")
        if isinstance(url, str) and url:
            return url
        return _extract_csv_url(data.get("data"))
    return ""


@safe.decorator(message="Failed to download CSV")
def _stream_csv_rows(
    client: "Client",
    url: str,
    rq: Optional[RangeQuery] = None,
    timeout: int = 30,
) -> Optional[List[Dict[str, Any]]]:
    """Stream-download CSV and parse rows, reusing the Client session (preserving proxy/retry config)."""
    import csv
    import time
    from collections import deque

    resp = client._session.get(url, stream=True, timeout=timeout)
    resp.raise_for_status()
    resp.encoding = "utf-8"

    lines = resp.iter_lines(decode_unicode=True)
    next(lines, None)  # skip header

    max_len = rq.tail if rq is not None and rq.tail is not None else None
    rows: Union[deque[Dict[str, Any]], List[Dict[str, Any]]] = deque(maxlen=max_len) if max_len is not None else []

    # Pre-compute timestamp lower bound for ``last`` mode
    last_start_ts: Optional[int] = None
    if rq is not None and rq.last is not None:
        last_start_ts = int(time.time() * 1000) - rq.last

    for row in csv.reader(lines):
        if not row or len(row) < 2:
            continue
        try:
            step = int(row[0])
            value = float(row[1])
        except (ValueError, IndexError):
            continue

        item: Dict[str, Any] = {"step": step, "value": value}

        if len(row) >= 3:
            try:
                item["timestamp"] = int(row[2])
            except (ValueError, IndexError):
                pass

        if rq is not None:
            # --- ``last`` mode: filter by timestamp >= (now - last) ---
            if last_start_ts is not None:
                ts = item.get("timestamp")
                if ts is None:
                    console.warning(f"CSV row missing timestamp column: {row}")
                    continue
                if ts < last_start_ts:
                    continue
            # --- timestamp range mode ---
            elif rq.type == "timestamp":
                ts = item.get("timestamp")
                if ts is None:
                    console.warning(f"CSV row missing timestamp column: {row}")
                    continue
                if rq.start is not None and ts < rq.start:
                    continue
                if rq.end is not None and ts > rq.end:
                    break
            # --- step range mode ---
            else:
                if rq.start is not None and step < rq.start:
                    continue
                if rq.end is not None and step > rq.end:
                    break
            # --- head early-stop ---
            if rq.head is not None and len(rows) >= rq.head:
                break

        rows.append(item)

    return list(rows)


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
        """指标数据列表。SCALAR 类型为采样后的折线点（含 index/data/timestamp），
        MEDIA 类型为 ``[{index, items}]`` 结构。

        .. note::
           当通过 ``Metrics`` 批量查询并使用 ``range_query`` 的 ``head`` / ``tail`` 时，
           ``head`` / ``tail`` 是 post-sampling 操作——在采样/下载完成后截取，而非对原始全量数据截取后再采样。
        """
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
        x_type: str = "step",
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> Dict[str, Any]:
        return {
            "projectId": project_id,
            "xType": x_type,
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

        # 1. 获取折线数据 — 使用 key-indexed lookup 保证对齐
        scalar_resp = self._post("/house/metrics/scalar", data=payload)
        if scalar_resp.ok and isinstance(scalar_resp.data, list):
            scalar_by_key = _align_entries_by_key(scalar_resp.data)
            res["metrics"] = scalar_by_key.get(self.key, {}).get("metrics", [])
        if not res.get("metrics"):
            return res

        # 2. 获取统计值 — step/time 并发，key-indexed 合并
        step_payload = {**payload, "xType": "step"}
        time_payload = {**payload, "xType": "timestamp"}
        step_resp, time_resp = self._concurrent_request(
            [
                (self._post, "/house/metrics/scalar/value", {"data": step_payload}),
                (self._post, "/house/metrics/scalar/value", {"data": time_payload}),
            ]
        )
        step_list = step_resp.data if step_resp.ok and isinstance(step_resp.data, list) else []
        time_list = time_resp.data if time_resp.ok and isinstance(time_resp.data, list) else []
        value_list = _merge_value_stats(step_list, time_list, [self.key])
        if value_list:
            for field in _SCALAR_STATISTIC_FIELDS:
                val = value_list[0].get(field)
                if val:
                    res[field] = val
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
    def _fetch_file_presigned_urls(entity: BaseEntity, paths: List[str]) -> Dict[str, str]:
        """通过完整资源路径批量获取预签名下载链接，返回 path → url 映射。"""
        if not paths:
            return {}
        resp = entity._post("/files/presigned/get", data={"paths": paths})
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

    def export_logs(self, start: int = 0, rows: int = 500_000) -> ApiResponseType:
        """导出日志数据为 .log 文件，返回对象存储中的文件信息。"""
        if self._metric_type != "LOG":
            return ApiResponseType(ok=False, errmsg="export_logs() only supports LOG metric_type", data=None)
        data: Dict[str, Any] = {
            "projectId": self._project_id,
            "experimentId": self._run_id,
            "start": start,
            "rows": min(rows, 500_000),
        }
        if self._root_pro_id:
            data["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            data["rootExpId"] = self._root_exp_id
        resp = self._post("/house/metrics/log/export", data=data)
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

    内部通过 ``_ensure_batch()`` 缓存结果，避免 ``__iter__`` 与 ``json()`` 重复请求。
    当 key 数量超过 ``_BATCH_SIZE``（默认 4）时自动分批，多批并发执行。

    .. note::
       ``range_query`` 的 ``head`` 和 ``tail`` 参数是 post-sampling 操作：
       在 sampled 模式下，服务端先做 LTTB 降采样，再对采样结果截取 head/tail；
       在 CSV 全量下载模式（``all`` 或 ``range_query``）下，先下载完整数据并做范围过滤，
       再对过滤后的结果截取 head/tail。因此 head/tail 不等同于对原始全量数据截取后再采样。

    用法::

        for m in experiment.metrics(keys=["loss", "acc"], metric_type="SCALAR"):
            print(m.key, m.metrics)
    """

    # 每批最大 key 数量，超出后自动拆分为多批并发
    _BATCH_SIZE = MAX_METRIC_KEY_BATCH_SIZE

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
        range_query: Optional[RangeQuery] = None,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> None:
        super().__init__(ctx)
        validate_metric_keys(keys)
        validate_metric_type(metric_type, keys[0])
        if metric_type == "LOG":
            raise ValueError("Metrics does not support LOG metric_type, use Experiment.logs() instead")
        if range_query is not None and metric_type != "SCALAR":
            raise ValueError("range_query is only supported for SCALAR metric_type")
        self._project_id = project_id
        self._run_id = run_id
        # 去重，保持插入顺序
        self._keys = list(dict.fromkeys(keys))
        self._metric_type = metric_type
        self._ignore_timestamp = ignore_timestamp
        self._media_step = media_step
        self._all = all
        self._range_query = range_query
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
        # 批量结果缓存：避免 __iter__ 与 json() 重复请求
        self._cached_list: Optional[List[Metric]] = None

    # ------------------------------------------------------------------
    # 公开接口
    # ------------------------------------------------------------------

    def _ensure_batch(self) -> List[Metric]:
        """Fetch (if needed) and cache batch Metric objects."""
        if self._cached_list is not None:
            return self._cached_list
        self._cached_list = list(self._fetch_batch())
        return self._cached_list

    def __iter__(self) -> Iterator[Metric]:
        yield from self._ensure_batch()

    def json(self) -> Dict[str, Any]:
        self._page_info["list"] = [
            {k: v for k, v in m.json().items() if k not in _METRIC_SHARED_KEYS} for m in self._ensure_batch()
        ]
        return self._page_info

    # ------------------------------------------------------------------
    # Batch dispatch: 分批限流
    # ------------------------------------------------------------------

    def _fetch_batch(self) -> Iterator[Metric]:
        """根据 metric_type 和模式分发到具体的获取方法。"""
        if self._metric_type == "SCALAR":
            if self._range_query is not None or self._all:
                data_list = self._batch_keys(self._fetch_scalar_csv)
            else:
                data_list = self._batch_keys(self._fetch_scalar_lines)
        else:
            # media 后端已支持 columns 批量，无需分批
            if self._all:
                data_list = self._fetch_media_all()
            else:
                data_list = self._fetch_media_data()

        for data in data_list:
            yield self._build_metric(data.get("key", ""), data)

    def _batch_keys(self, fetch_fn: Callable[[List[str]], List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        """将 keys 按 ``_BATCH_SIZE`` 分批；单批直接执行，多批并发。"""
        keys = self._keys
        if len(keys) <= self._BATCH_SIZE:
            return fetch_fn(keys)

        chunks = [keys[i : i + self._BATCH_SIZE] for i in range(0, len(keys), self._BATCH_SIZE)]
        with SafeThreadPoolExecutor(max_workers=min(len(chunks), self._BATCH_SIZE)) as pool:
            futures = [pool.submit(fetch_fn, chunk) for chunk in chunks]
            results: List[Dict[str, Any]] = []
            for f in futures:
                results.extend(f.result())
            return results

    # ------------------------------------------------------------------
    # Metric 对象构建
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # 共享辅助方法
    # ------------------------------------------------------------------

    def _empty_scalar_results(self, keys: List[str]) -> List[Dict[str, Any]]:
        """返回 per-key 空结果列表（用于 early return）。"""
        return [{"projectId": self._project_id, "experimentId": self._run_id, "key": k, "metrics": []} for k in keys]

    def _build_value_stats_requests(self, keys: List[str]) -> List[tuple]:
        """构建 step/time value stats 的并发请求列表（2 路并发）。"""
        value_path = "/house/metrics/scalar/value"
        return [
            (
                self._post,
                value_path,
                {
                    "data": Metric._build_scalar_payload(
                        self._project_id,
                        self._run_id,
                        keys,
                        self._sample,
                        x_type=x_type,
                        root_pro_id=self._root_pro_id,
                        root_exp_id=self._root_exp_id,
                    )
                },
            )
            for x_type in ("step", "timestamp")
        ]

    @staticmethod
    def _extract_value_stats(
        step_resp: ApiResponseType,
        time_resp: ApiResponseType,
        keys: List[str],
    ) -> List[Dict[str, Any]]:
        """从并发的 step/time 响应中提取并合并 value stats。"""
        step_list = step_resp.data if step_resp.ok and isinstance(step_resp.data, list) else []
        time_list = time_resp.data if time_resp.ok and isinstance(time_resp.data, list) else []
        return _merge_value_stats(step_list, time_list, keys)

    # ------------------------------------------------------------------
    # Scalar: 折线数据 + 统计值 (后端 columns 批量)
    # ------------------------------------------------------------------

    def _fetch_scalar_lines(self, keys: List[str]) -> List[Dict[str, Any]]:
        """获取标量折线数据 + step/time 统计值，3 路并发。

        后端 ``POST /house/metrics/scalar`` 和 ``/scalar/value`` 的 ``columns``
        数组天然支持多 key，此处将 keys 打包为一个批量请求。
        """
        # 3 路并发：折线数据 + step 统计 + time 统计
        requests: List[tuple] = [
            (
                self._post,
                "/house/metrics/scalar",
                {
                    "data": Metric._build_scalar_payload(
                        self._project_id,
                        self._run_id,
                        keys,
                        self._sample,
                        root_pro_id=self._root_pro_id,
                        root_exp_id=self._root_exp_id,
                    )
                },
            ),
        ]
        requests.extend(self._build_value_stats_requests(keys))

        scalar_resp, step_resp, time_resp = self._concurrent_request(requests)

        scalar_list = scalar_resp.data if scalar_resp.ok and isinstance(scalar_resp.data, list) else []
        scalar_by_key = _align_entries_by_key(scalar_list)
        metrics_by_key: Dict[str, Any] = {key: scalar_by_key.get(key, {}).get("metrics", []) for key in keys}
        value_list = self._extract_value_stats(step_resp, time_resp, keys)
        value_by_key: Dict[str, Dict[str, Any]] = {keys[i]: v for i, v in enumerate(value_list)}

        return self._build_scalar_results(keys, metrics_by_key, value_by_key)

    # ------------------------------------------------------------------
    # Scalar: CSV 全量下载 + 统计值 (range_query 或 all 模式)
    # ------------------------------------------------------------------

    def _fetch_scalar_csv(self, keys: List[str]) -> List[Dict[str, Any]]:
        """CSV 全量下载 + value stats，URL 获取与下载均并发。

        value stats 通过后端 ``columns`` 批量获取；
        CSV presigned URL 通过 ``GET /experiment/{run_id}/column/csv`` per-key 获取，
        多 key 时并发拉取 URL + 并发下载 CSV。
        """
        # 并发：step 统计 + time 统计 + per-key CSV URL
        requests: List[tuple] = self._build_value_stats_requests(keys)
        for key in keys:
            requests.append((self._get, f"/experiment/{self._run_id}/column/csv", {"params": {"key": key}}))

        all_resps = self._concurrent_request(requests)
        step_resp, time_resp = all_resps[0], all_resps[1]
        csv_resps = all_resps[2:]

        value_list = self._extract_value_stats(step_resp, time_resp, keys)
        value_by_key: Dict[str, Dict[str, Any]] = {keys[i]: v for i, v in enumerate(value_list)}

        # 提取 presigned CSV URL（按 key 索引）
        url_by_key: Dict[str, Optional[str]] = {}
        for i, key in enumerate(keys):
            url = ""
            if i < len(csv_resps) and csv_resps[i].ok and csv_resps[i].data:
                url = _extract_csv_url(csv_resps[i].data)
            url_by_key[key] = url or None

        # 并发下载 CSV 并解析
        urls_ordered = [url_by_key.get(key) for key in keys]
        csv_rows_list = self._download_csvs(urls_ordered)
        metrics_by_key: Dict[str, Any] = {
            keys[i]: csv_rows_list[i] if i < len(csv_rows_list) else [] for i, key in enumerate(keys)
        }

        return self._build_scalar_results(keys, metrics_by_key, value_by_key)

    # ------------------------------------------------------------------
    # Scalar results 构建
    # ------------------------------------------------------------------

    def _build_scalar_results(
        self,
        keys: List[str],
        metrics_by_key: Dict[str, Any],
        value_by_key: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """将折线/CSV 数据与 value stats 合并为 per-key 结果字典。

        所有数据通过 key 索引查找，保证与请求 keys 顺序对齐，不依赖后端返回顺序。
        """
        results: List[Dict[str, Any]] = []
        for key in keys:
            data: Dict[str, Any] = {
                "projectId": self._project_id,
                "experimentId": self._run_id,
                "key": key,
                "metrics": metrics_by_key.get(key, []),
            }
            stats = value_by_key.get(key, {})
            if stats:
                data.update(stats)
            results.append(data)
        return results

    # ------------------------------------------------------------------
    # CSV 并发下载
    # ------------------------------------------------------------------

    def _download_csvs(self, urls: List[Optional[str]]) -> List[List[Dict[str, Any]]]:
        """通过线程池并发下载并解析 CSV 文件。"""
        if not any(urls):
            return [[] for _ in urls]

        with SafeThreadPoolExecutor(max_workers=min(len(urls), self._BATCH_SIZE)) as pool:
            futures = [
                pool.submit(_stream_csv_rows, self._ctx.client, url, self._range_query) if url else None for url in urls
            ]
            return [(f.result() or []) if f is not None else [] for f in futures]

    # ------------------------------------------------------------------
    # Media fetch（后端 columns 批量，单次请求，无需分批）
    # ------------------------------------------------------------------

    def _fetch_media_data(self) -> List[Dict[str, Any]]:
        """获取媒体数据（单步），后端 columns 批量一次返回。"""
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
            return self._empty_scalar_results(self._keys)
        resp_data = raw_resp.data
        if not isinstance(resp_data, dict):
            return self._empty_scalar_results(self._keys)

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
        results: List[Dict[str, Any]] = []
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
            results.append(data)
        return results

    def _fetch_media_all(self) -> List[Dict[str, Any]]:
        """获取全部媒体数据，后端 columns 批量一次返回。"""
        payload = Metric._build_media_payload(
            self._project_id,
            self._run_id,
            self._keys,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/f_media", data=payload)
        if not raw_resp.ok or not raw_resp.data:
            return self._empty_scalar_results(self._keys)
        raw_list = raw_resp.data
        if not isinstance(raw_list, list):
            return self._empty_scalar_results(self._keys)

        prefix = f"{self._project_id}/{self._run_id}"
        all_paths = [p for entry in raw_list for m in entry.get("metrics", []) for p in m.get("data", [])]
        url_map = Metric._fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched (all): run_id[{self._run_id}] - {len(all_paths)} items across {len(self._keys)} keys, requesting presigned urls..."
            )

        key_to_entry: Dict[str, Dict[str, Any]] = {e.get("key", ""): e for e in raw_list}
        results: List[Dict[str, Any]] = []
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
            results.append(data)
        return results
