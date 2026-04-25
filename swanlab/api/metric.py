"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/20
@description: Metric 实体类 — 指标序列的查询与操作
"""

from typing import Any, Dict, List, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiColumnCsvExportType, ApiResponseType
from swanlab.api.typings.metric import ApiLogSeriesType, ApiMediaSeriesType, ApiMediaType, ApiScalarSeriesType
from swanlab.api.utils import get_properties, validate_metric_log_level, validate_metric_type


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
        log_level: str = "INFO",
        metric_type: str = "SCALAR",
        data: Optional[Dict[str, Any]] = None,
        ignore_timestamp: bool = False,
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
        # TODO: 采样值， scalar 时生效，logs 时降级到 1000
        self._sample = sample
        # 偏移量，仅对 Log metric_type 有效， 默认为 0
        self._offset = log_offset
        self._log_level = log_level

    # 类型 → 加载方法 的分发表，新增类型只需在此注册
    _FETCH_DISPATCH = {
        "SCALAR": "_fetch_scalar",
        "MEDIA": "_fetch_media",
        "LOG": "_fetch_logs",
    }

    def _ensure_data(self) -> Dict[str, Any]:
        if self._data is None:
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

    # ------------------------------------------------------------------
    # 请求辅助函数
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_first(resp: ApiResponseType) -> Optional[Dict[str, Any]]:
        """从列表型 API 响应中提取第一个元素，失败返回 None。"""
        if resp.ok and isinstance(resp.data, list) and resp.data:
            return resp.data[0]
        return None

    def _build_scalar_payload(self) -> Dict[str, Any]:
        return {
            "projectId": self.project_id,
            "xType": "step",
            "range": [0, 0],
            "columns": [{"experimentId": self.run_id, "key": self.key}],
        }

    def _build_media_payload(self) -> Dict[str, Any]:
        return {
            "projectId": self.project_id,
            "columns": [{"experimentId": self.run_id, "key": self.key}],
        }

    def _build_log_params(self) -> Dict[str, Any]:
        return {
            "projectId": self.project_id,
            "experimentId": self.run_id,
            "size": 1000,  # 硬编码为 1000
            "epoch": self._offset,
            "level": self._log_level,
        }

    # ------------------------------------------------------------------
    # 类型专属加载
    # ------------------------------------------------------------------

    def _fetch_scalar(self) -> ApiScalarSeriesType:
        res = ApiScalarSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = self._build_scalar_payload()

        # 1. 获取折线数据
        raw_data = self._extract_first(self._post("/house/metrics/scalar", data=payload))
        if raw_data is None:
            return res
        res["metrics"] = raw_data.get("metrics", [])

        # 2. 获取统计值
        stat_data = self._extract_first(self._post("/house/metrics/scalar/value", data=payload))
        if stat_data is None:
            return res
        for field in ("min", "max", "avg", "median", "latest"):
            res[field] = stat_data.get(field, {})
        return res

    def _fetch_media(self) -> ApiMediaSeriesType:
        res = ApiMediaSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        payload = self._build_media_payload()
        raw_resp = self._post("/house/metrics/f_media", data=payload)
        raw_data = self._extract_first(raw_resp)
        if raw_data is None:
            return res
        metrics: List[ApiMediaType] = []
        prefix = f"{self.project_id}/{self.run_id}"
        for entry in raw_data.get("metrics", []):
            paths = entry.get("data", [])
            mores = entry.get("more", [])
            items = []
            for i, path in enumerate(paths):
                item = {"path": path}
                if i < len(mores) and isinstance(mores[i], dict):
                    item.update(mores[i])
                items.append(item)
            metrics.append({"index": entry.get("index", 0), "prefix": prefix, "items": items})

        res["metrics"] = metrics
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
        """
        导出列数据为 CSV。

        :return: ApiResponseType，成功时 data 包含临时下载 URL
        """
        if self.metric_type != "SCALAR":
            err_msg = "export_csv() only support SCALAR metric_type"
            return ApiResponseType(ok=False, errmsg=err_msg, data=None)
        resp = self._get(f"/experiment/{self._run_id}/column/csv", params={"key": self.key})
        if not resp.ok:
            return resp

        data = resp.data
        if isinstance(data, list) and data:
            url = data[0].get("url", "")
        elif isinstance(data, dict):
            url = data.get("url", "")
        else:
            return ApiResponseType(ok=False, errmsg="Invalid response format", data=None)
        return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))

    def json(self) -> Dict[str, Any]:
        result = get_properties(self)
        data = self._ensure_data()

        if self._metric_type == "SCALAR":
            for field in ("min", "max", "avg", "median", "latest"):
                val = data.get(field)
                if val:
                    result[field] = val

        if self._metric_type == "LOG":
            result.pop("metrics", None)
        else:
            result.pop("logs", None)
            result.pop("count", None)

        if self._ignore_timestamp:
            timestamp_items = result.get("metrics", []) or result.get("logs", [])
            for item in cast(List[Dict[str, Any]], timestamp_items):
                if isinstance(item, dict):
                    item.pop("timestamp", None)

        return result
