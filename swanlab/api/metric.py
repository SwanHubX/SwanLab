"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/20
@description: Metric 实体类 — 指标序列的查询与操作
"""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.helper import (
    SCALAR_STATISTIC_FIELDS,
    STEP_AXIS,
    align_entries_by_key,
    axis_request_params,
    build_export_payload,
    build_media_items,
    build_media_payload,
    build_scalar_payload,
    builtin_x_axis,
    extract_first,
    extract_value_stats,
    fetch_file_presigned_urls,
    fetch_presigned_urls,
    get_properties,
    group_keys_by_axis,
    merge_value_stats,
    stream_export_csv,
    validate_metric_keys,
    validate_metric_log_level,
    validate_metric_type,
    validate_x_axis,
)
from swanlab.api.typings import ApiColumnCsvExportType, ApiResponseType
from swanlab.api.typings.common import (
    MAX_CONCURRENT_COUNT,
    MAX_METRIC_KEY_BATCH_SIZE,
    ApiMetricKeyTypeLiteral,
    ApiMetricLogLevelLiteral,
    RangeQuery,
)
from swanlab.api.typings.metric import (
    ApiLogSeriesType,
    ApiMediaSeriesType,
    ApiMetricXAxisType,
    ApiScalarSeriesType,
)
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor

_METRIC_SHARED_KEYS = frozenset({"project_id", "run_id", "metric_type"})


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
        created_at: int,
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
        experiment_name: str = "",
        x_axis: str = "auto",
        proj_path: str = "",
    ) -> None:
        super().__init__(ctx)
        validate_metric_type(metric_type, key)
        if metric_type == "LOG":
            validate_metric_log_level(log_level)
        validate_x_axis(x_axis, metric_type=metric_type)
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
        self._created_at = created_at
        self._experiment_name = experiment_name
        # X 轴：auto 按 DEFAULT 视图配置逐 key 解析（需 proj_path，缺失时回落 step）
        self._x_axis = x_axis
        self._proj_path = proj_path

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
        # createdAt 为 House 查询的数据入库时间下界，必传
        params["createdAt"] = self._created_at
        return params

    # ------------------------------------------------------------------
    # 类型专属加载
    # ------------------------------------------------------------------

    def _resolve_axis(self) -> ApiMetricXAxisType:
        """解析当前 key 的 X 轴：显式字面量直接映射，auto 走 DEFAULT 视图发现（缺 proj_path 回落 step）。"""
        builtin = builtin_x_axis(self._x_axis)
        if builtin is not None:
            return builtin
        if self._x_axis == "auto":
            if self._proj_path and self._metric_type == "SCALAR":
                from swanlab.api.view import fetch_default_view_xaxis_map

                axis_map = fetch_default_view_xaxis_map(self._ctx, self._proj_path)
                return axis_map.get(self.key) or STEP_AXIS
            return STEP_AXIS
        return {"type": "CUSTOM", "key": self._x_axis}

    def _fetch_scalar(self) -> ApiScalarSeriesType:
        res = ApiScalarSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        if self._root_pro_id:
            res["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            res["rootExpId"] = self._root_exp_id

        axis = self._resolve_axis()
        res["xAxis"] = axis
        x_type, x_key = axis_request_params(axis)
        payload = build_scalar_payload(
            self.project_id,
            self.run_id,
            self._created_at,
            [self.key],
            self._sample,
            x_type=x_type,
            x_key=x_key or None,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )

        # 1. 获取折线数据 — 使用 key-indexed lookup 保证对齐
        scalar_resp = self._post("/house/metrics/scalar", data=payload)
        if scalar_resp.ok and isinstance(scalar_resp.data, list):
            scalar_by_key = align_entries_by_key(scalar_resp.data)
            res["metrics"] = scalar_by_key.get(self.key, {}).get("metrics", [])
        if not res.get("metrics"):
            return res

        # 2. 获取统计值 — step/time 并发，key-indexed 合并
        #    /scalar/value 不支持 xKey，stats payload 一律新建（不能从折线 payload 派生）
        step_payload = build_scalar_payload(
            self.project_id,
            self.run_id,
            self._created_at,
            [self.key],
            self._sample,
            x_type="step",
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        time_payload = build_scalar_payload(
            self.project_id,
            self.run_id,
            self._created_at,
            [self.key],
            self._sample,
            x_type="timestamp",
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        step_resp, time_resp = self._concurrent_request(
            [
                (self._post, "/house/metrics/scalar/value", {"data": step_payload}),
                (self._post, "/house/metrics/scalar/value", {"data": time_payload}),
            ]
        )
        step_list = step_resp.data if step_resp.ok and isinstance(step_resp.data, list) else []
        time_list = time_resp.data if time_resp.ok and isinstance(time_resp.data, list) else []
        value_list = merge_value_stats(step_list, time_list, [self.key])
        if value_list:
            for field in SCALAR_STATISTIC_FIELDS:
                val = value_list[0].get(field)
                if val:
                    res[field] = val
        return res

    def _fetch_media(self) -> ApiMediaSeriesType:
        res = ApiMediaSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = build_media_payload(
            self.project_id,
            self.run_id,
            self._created_at,
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

        res["steps"] = data.get("steps") or []
        step_val = data.get("step")
        if step_val is not None:
            res["step"] = step_val

        # metrics 可能为 null（该步无数据）或含 None 条目，过滤后再匹配 key
        metrics_raw: List[Dict[str, Any]] = [m for m in (data.get("metrics") or []) if isinstance(m, dict)]
        metric_entry = next((m for m in metrics_raw if m.get("key") == self.key), None)
        if metric_entry is None:
            return res

        prefix = f"{self.project_id}/{self.run_id}"
        all_paths = metric_entry.get("data") or []
        url_map = fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched: run_id[{self.run_id}], key[{self.key}] - {len(all_paths)} items, requesting presigned urls..."
            )
        items = build_media_items(metric_entry, url_map)
        res["metrics"] = [{"index": data.get("step") or 0, "items": items}]
        return res

    def _fetch_media_all(self) -> ApiMediaSeriesType:
        res = ApiMediaSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = build_media_payload(
            self.project_id,
            self.run_id,
            self._created_at,
            [self.key],
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        raw_resp = self._post("/house/metrics/f_media", data=payload)
        raw_data = extract_first(raw_resp)
        if raw_data is None:
            return res

        prefix = f"{self.project_id}/{self.run_id}"
        # metrics 可能为 null 或含 None 条目，过滤后再展开
        entries = [e for e in (raw_data.get("metrics") or []) if isinstance(e, dict)]
        all_paths = [p for entry in entries for p in (entry.get("data") or [])]
        url_map = fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched (all): run_id[{self.run_id}], key[{self.key}] - {len(all_paths)} items, requesting presigned urls..."
            )
        res["metrics"] = [
            {"index": entry.get("index", 0), "items": build_media_items(entry, url_map)} for entry in entries
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
        """导出列数据为 CSV。

        通过 ``POST /house/metrics/scalar/export`` 获取 ``key``，
        再通过 ``/files/presigned/get`` 转换为预签名下载链接。
        """
        if self.metric_type != "SCALAR":
            return ApiResponseType(ok=False, errmsg="export_csv() only support SCALAR metric_type", data=None)
        payload = build_export_payload(
            self._project_id,
            self._run_id,
            self._created_at,
            [self.key],
            experiment_name=self._experiment_name,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        resp = self._post("/house/metrics/scalar/export", data=payload)
        if not resp.ok or not resp.data:
            return resp
        cos_key = resp.data.get("cosKey", "") if isinstance(resp.data, dict) else ""
        if not cos_key:
            return ApiResponseType(ok=False, errmsg="Invalid response format: missing cosKey", data=None)
        url_map = fetch_file_presigned_urls(self, [cos_key])
        url = url_map.get(cos_key, "")
        if not url:
            return ApiResponseType(ok=False, errmsg="Failed to get presigned download URL", data=None)
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
        # createdAt 为 House 查询的数据入库时间下界，必传
        data["createdAt"] = self._created_at
        resp = self._post("/house/metrics/log/export", data=data)
        if not resp.ok or not resp.data:
            return resp
        cos_key = resp.data.get("cosKey", "")
        if not cos_key:
            return ApiResponseType(ok=False, errmsg="Invalid response format: missing cosKey", data=None)
        url_map = fetch_file_presigned_urls(self, [cos_key])
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
            # 回显 metrics[].index 的实际轴语义（自定义 x / 内置轴 / 回落 step）
            if "xAxis" in data:
                result["xAxis"] = data["xAxis"]
            for field in SCALAR_STATISTIC_FIELDS:
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
        metric_type: ApiMetricKeyTypeLiteral,
        sample: int = 1500,
        ignore_timestamp: bool = False,
        media_step: Optional[int] = None,
        all: bool = False,
        range_query: Optional[RangeQuery] = None,
        root_pro_id: str = "",
        root_exp_id: str = "",
        created_at: int,
        experiment_name: str = "",
        x_axis: str = "auto",
        proj_path: str = "",
    ) -> None:
        super().__init__(ctx)
        validate_metric_keys(keys)
        validate_metric_type(metric_type, keys[0])
        if metric_type == "LOG":
            raise ValueError("Metrics does not support LOG metric_type, use Experiment.logs() instead")
        if range_query is not None and metric_type != "SCALAR":
            raise ValueError("range_query is only supported for SCALAR metric_type")
        validate_x_axis(x_axis, metric_type=metric_type)
        # CSV 全量管线围绕自定义 x 列构建：显式内置时间轴 + all/range_query 在构造期报错，
        # 不静默回落 step（回落会掩盖"用户要的轴没生效"）；auto 的 relative_time 分组在
        # 分组解析期回落 step + warning（见 _fetch_scalar_csv_grouped）
        if (range_query is not None or all) and x_axis in ("time", "relative_time"):
            raise ValueError(
                f"x_axis={x_axis!r} is not supported in CSV mode (all=True or range_query); "
                "use the sampled mode (default) or a custom x column key instead"
            )
        # type="custom" 的构造期校验：显式内置轴时所有分组必然非 custom，直接报错；
        # auto 模式需等发现 + 分组解析后按分组判定（见 _validate_custom_range_groups）
        if range_query is not None and range_query.type == "custom":
            builtin = builtin_x_axis(x_axis)
            if builtin is not None and builtin.get("type") != "CUSTOM":
                raise ValueError(
                    "range_query type='custom' filters on custom x-axis values and requires a custom x axis; "
                    f"got x_axis={x_axis!r}. Pass x_axis=<custom column key> instead"
                )
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
        self._created_at = created_at
        self._experiment_name = experiment_name
        self._x_axis = x_axis
        self._proj_path = proj_path
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
                data_list = self._fetch_scalar_csv_grouped()
            else:
                data_list = self._fetch_scalar_sampled_grouped()
        else:
            # media 后端已支持 columns 批量，无需分批
            if self._all:
                data_list = self._fetch_media_all()
            else:
                data_list = self._fetch_media_data()

        for data in data_list:
            yield self._build_metric(data.get("key", ""), data)

    # ------------------------------------------------------------------
    # X 轴解析与分组（单次 House 请求只能携带一个轴）
    # ------------------------------------------------------------------

    def _resolve_axes(self) -> Dict[str, ApiMetricXAxisType]:
        """逐 key 解析 X 轴：显式字面量统一映射；auto 走 DEFAULT 视图发现，未知 key 回落 step。"""
        builtin = builtin_x_axis(self._x_axis)
        if builtin is not None:
            return {key: builtin for key in self._keys}
        if self._x_axis == "auto":
            if self._proj_path:
                from swanlab.api.view import fetch_default_view_xaxis_map

                axis_map = fetch_default_view_xaxis_map(self._ctx, self._proj_path)
                return {key: axis_map.get(key) or STEP_AXIS for key in self._keys}
            return {key: STEP_AXIS for key in self._keys}
        custom: ApiMetricXAxisType = {"type": "CUSTOM", "key": self._x_axis}
        return {key: custom for key in self._keys}

    def _validate_custom_range_groups(self, axes: Dict[str, ApiMetricXAxisType]) -> None:
        """``type="custom"`` 合法性按分组判定（auto 模式的分组解析期校验）。

        任何回落 step（或解析为 SYSTEM 内置轴）的分组都不能按 x 值域过滤——静默按 step
        过滤即是"静默重解释"，直接报错并列出违规 keys。
        """
        if self._range_query is None or self._range_query.type != "custom":
            return
        violating = [key for key in self._keys if axes[key].get("type") != "CUSTOM"]
        if violating:
            raise ValueError(
                f"range_query type='custom' filters on custom x-axis values, but keys {violating} did not "
                "resolve to a custom x axis; pass an explicit x_axis=<custom column key> or fix the "
                "DEFAULT view chart config"
            )

    def _merge_grouped_results(self, merged: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
        """按调用方 key 顺序合并各分组结果。"""
        return [merged[key] for key in self._keys if key in merged]

    def _fetch_scalar_sampled_grouped(self) -> List[Dict[str, Any]]:
        """采样路径：按轴分组，组内走现有分批；结果按调用方 key 顺序合并。"""
        axes = self._resolve_axes()
        self._validate_custom_range_groups(axes)
        merged: Dict[str, Dict[str, Any]] = {}
        for axis, keys in group_keys_by_axis(axes):
            x_type, x_key = axis_request_params(axis)
            for data in self._batch_keys(keys, x_type=x_type, x_key=x_key):
                merged[data["key"]] = {**data, "xAxis": axis}
        return self._merge_grouped_results(merged)

    def _fetch_scalar_csv_grouped(self) -> List[Dict[str, Any]]:
        """CSV 全量路径：auto 解析出 SYSTEM relative_time 的分组回落 step + warning（v4）。

        显式 time/relative_time 已在构造期报错；回落信息无损（每个 item 仍携带 timestamp，
        relative time 可自行推导），xAxis 回显如实标注 step。
        """
        axes = self._resolve_axes()
        csv_axes: Dict[str, ApiMetricXAxisType] = {}
        for key, axis in axes.items():
            if axis.get("type") == "SYSTEM" and axis.get("key") == "relative_time":
                console.warning(
                    f"Metric key '{key}' is configured with a relative_time x-axis in the DEFAULT view, "
                    "which the CSV mode (all/range_query) does not support; falling back to 'step'. "
                    "Each point still carries its 'timestamp' for deriving relative time."
                )
                axis = STEP_AXIS
            csv_axes[key] = axis
        self._validate_custom_range_groups(csv_axes)
        merged: Dict[str, Dict[str, Any]] = {}
        for axis, keys in group_keys_by_axis(csv_axes):
            x_key = axis.get("key", "") if axis.get("type") == "CUSTOM" else ""
            for data in self._fetch_scalar_csv(keys, x_key=x_key):
                merged[data["key"]] = {**data, "xAxis": axis}
        return self._merge_grouped_results(merged)

    def _batch_keys(
        self,
        keys: List[str],
        x_type: str = "step",
        x_key: str = "",
    ) -> List[Dict[str, Any]]:
        """将 keys 按 ``_BATCH_SIZE`` 分批；单批直接执行，多批并发。"""
        if len(keys) <= self._BATCH_SIZE:
            return self._fetch_scalar_lines(keys, x_type=x_type, x_key=x_key)

        chunks = [keys[i : i + self._BATCH_SIZE] for i in range(0, len(keys), self._BATCH_SIZE)]
        with SafeThreadPoolExecutor(max_workers=min(len(chunks), self._BATCH_SIZE)) as pool:
            futures = [pool.submit(self._fetch_scalar_lines, chunk, x_type, x_key) for chunk in chunks]
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
            created_at=self._created_at,
            experiment_name=self._experiment_name,
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
                    "data": build_scalar_payload(
                        self._project_id,
                        self._run_id,
                        self._created_at,
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

    # ------------------------------------------------------------------
    # Scalar: 折线数据 + 统计值 (后端 columns 批量)
    # ------------------------------------------------------------------

    def _fetch_scalar_lines(self, keys: List[str], x_type: str = "step", x_key: str = "") -> List[Dict[str, Any]]:
        """获取标量折线数据 + step/time 统计值，3 路并发。

        后端 ``POST /house/metrics/scalar`` 和 ``/scalar/value`` 的 ``columns``
        数组天然支持多 key，此处将 keys 打包为一个批量请求。
        同一分组内的 keys 共享一个轴（xKey / xType）；stats 请求不支持 xKey，不带。
        """
        # 3 路并发：折线数据 + step 统计 + time 统计
        requests: List[tuple] = [
            (
                self._post,
                "/house/metrics/scalar",
                {
                    "data": build_scalar_payload(
                        self._project_id,
                        self._run_id,
                        self._created_at,
                        keys,
                        self._sample,
                        x_type=x_type,
                        x_key=x_key or None,
                        root_pro_id=self._root_pro_id,
                        root_exp_id=self._root_exp_id,
                    )
                },
            ),
        ]
        requests.extend(self._build_value_stats_requests(keys))

        scalar_resp, step_resp, time_resp = self._concurrent_request(requests)

        scalar_list = scalar_resp.data if scalar_resp.ok and isinstance(scalar_resp.data, list) else []
        scalar_by_key = align_entries_by_key(scalar_list)
        metrics_by_key: Dict[str, Any] = {key: scalar_by_key.get(key, {}).get("metrics", []) for key in keys}
        value_list = extract_value_stats(step_resp, time_resp, keys)
        value_by_key: Dict[str, Dict[str, Any]] = {keys[i]: v for i, v in enumerate(value_list)}

        return self._build_scalar_results(keys, metrics_by_key, value_by_key)

    # ------------------------------------------------------------------
    # Scalar: CSV 全量下载 + 统计值 (range_query 或 all 模式)
    # ------------------------------------------------------------------

    def _fetch_scalar_csv(self, keys: List[str], x_key: str = "") -> List[Dict[str, Any]]:
        """CSV 全量下载 + value stats，使用 House 批量导出接口。

        将 keys 按 ``MAX_METRIC_KEY_BATCH_SIZE``（16）个一组分批，4 线程并发获取。
        自定义 x 轴下 ``x_key`` 作为额外的普通导出列追加，占用一个分批槽位
        （每批 y keys 上限降为 15）。
        通过 ``POST /house/metrics/scalar/export`` 一次性导出每批 key 到 CSV 文件，
        再通过 ``/files/presigned/get`` 获取预签名下载链接。
        value stats 批量获取，与导出请求并发执行。
        """
        chunk_size = MAX_METRIC_KEY_BATCH_SIZE - 1 if x_key else MAX_METRIC_KEY_BATCH_SIZE
        if len(keys) <= chunk_size:
            return self._fetch_scalar_csv_batch(keys, x_key=x_key)

        chunks = [keys[i : i + chunk_size] for i in range(0, len(keys), chunk_size)]
        with SafeThreadPoolExecutor(max_workers=MAX_CONCURRENT_COUNT) as pool:
            futures = [pool.submit(self._fetch_scalar_csv_batch, chunk, x_key) for chunk in chunks]
            results: List[Dict[str, Any]] = []
            for f in futures:
                results.extend(f.result())
            return results

    def _fetch_scalar_csv_batch(self, keys: List[str], x_key: str = "") -> List[Dict[str, Any]]:
        """单批 CSV 导出 + value stats（≤16 列；自定义 x 时 ≤15 y keys + 1 x 列）。

        通过 ``POST /house/metrics/scalar/export`` 一次性导出所有 key 到一个 CSV 文件，
        再通过 ``/files/presigned/get`` 获取预签名下载链接。
        value stats 批量获取，与导出请求并发执行；stats 请求只含 y keys，不含 x 列。
        """
        # 并发：step 统计 + time 统计 + 批量 CSV 导出
        requests: List[tuple] = self._build_value_stats_requests(keys)
        export_columns = keys + ([x_key] if x_key else [])
        export_payload = build_export_payload(
            self._project_id,
            self._run_id,
            self._created_at,
            export_columns,
            experiment_name=self._experiment_name,
            root_pro_id=self._root_pro_id,
            root_exp_id=self._root_exp_id,
        )
        requests.append((self._post, "/house/metrics/scalar/export", {"data": export_payload}))

        all_resps = self._concurrent_request(requests)
        step_resp, time_resp = all_resps[0], all_resps[1]
        export_resp = all_resps[2]

        value_list = extract_value_stats(step_resp, time_resp, keys)
        value_by_key: Dict[str, Dict[str, Any]] = {keys[i]: v for i, v in enumerate(value_list)}

        # cosKey → presigned URL → download → parse per-key rows
        metrics_by_key: Dict[str, Any] = {key: [] for key in keys}
        cos_key = ""
        if export_resp.ok and isinstance(export_resp.data, dict):
            cos_key = export_resp.data.get("cosKey", "")
        if cos_key:
            url_map = fetch_file_presigned_urls(self, [cos_key])
            url = url_map.get(cos_key, "")
            if url:
                parsed = stream_export_csv(self._ctx.client, url, keys, self._range_query, x_key=x_key)
                if parsed:
                    metrics_by_key = parsed

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
    # Media fetch（后端 columns 批量，单次请求，无需分批）
    # ------------------------------------------------------------------

    def _fetch_media_data(self) -> List[Dict[str, Any]]:
        """获取媒体数据（单步），后端 columns 批量一次返回。"""
        payload = build_media_payload(
            self._project_id,
            self._run_id,
            self._created_at,
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

        steps = resp_data.get("steps") or []
        current_step = resp_data.get("step")
        # metrics 可能为 null（请求的步数无数据）或含 None 条目，过滤后再展开
        metrics_raw: List[Dict[str, Any]] = [m for m in (resp_data.get("metrics") or []) if isinstance(m, dict)]

        prefix = f"{self._project_id}/{self._run_id}"
        all_paths = [p for entry in metrics_raw for p in (entry.get("data") or [])]
        url_map = fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
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
                items = build_media_items(entry, url_map)
                data["metrics"] = [{"index": current_step or 0, "items": items}]
            results.append(data)
        return results

    def _fetch_media_all(self) -> List[Dict[str, Any]]:
        """获取全部媒体数据，后端 columns 批量一次返回。"""
        payload = build_media_payload(
            self._project_id,
            self._run_id,
            self._created_at,
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
        # 列表元素可能为 None（无数据的 key），metrics 字段也可能为 null，逐层过滤
        entries = [e for e in raw_list if isinstance(e, dict)]
        all_paths = [
            p
            for entry in entries
            for m in (entry.get("metrics") or [])
            if isinstance(m, dict)
            for p in (m.get("data") or [])
        ]
        url_map = fetch_presigned_urls(self, prefix, all_paths) if all_paths else {}
        if all_paths:
            console.debug(
                f"Media fetched (all): run_id[{self._run_id}] - {len(all_paths)} items across {len(self._keys)} keys, requesting presigned urls..."
            )

        key_to_entry: Dict[str, Dict[str, Any]] = {e.get("key", ""): e for e in entries}
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
                for m in entry.get("metrics") or []:
                    if not isinstance(m, dict):
                        continue
                    items = build_media_items(m, url_map)
                    metrics_list.append({"index": m.get("index", 0), "items": items})
                data["metrics"] = metrics_list
            results.append(data)
        return results
