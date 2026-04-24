"""
@author: caddiesnew
@file: column.py
@time: 2026/4/20
@description: Column 实体类 — 实验列的查询与操作
"""

from typing import Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiColumnCsvExportType, ApiMetricTypeLiteral, ApiResponseType
from swanlab.api.typings.metric import ApiLogType, ApiMediaType, ApiMetricType, ApiScalarSeriesType, ApiScalarType
from swanlab.api.utils import get_properties, resovle_run_path, validate_column_params, validate_metric_type


class Metric(BaseEntity):
    """
    表示一个 SwanLab 指标列 (非单个数值，而是一组序列)
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        project_id: str,
        run_id: str,
        key: Optional[str] = "",
        sample: int = 1500,
        metric_type: str = "SCALAR",
        data: Optional[Any] = None,
    ) -> None:
        super().__init__(ctx)
        validate_metric_type(metric_type, key)
        self._project_id = project_id
        self._run_id = run_id
        self._key = key
        self._data = data
        self._metric_type = metric_type

        # TODO: 采样值，仅在 scalar 时生效， 待接入
        self._sample = sample

    def _ensure_data(self) -> Dict[str, Any]:
        if self._data is None:
            if self._metric_type == "SCALAR":
                self._data = self._fetch_scalar()
                print(self._data)
            elif self._metric_type == "MEDIA":
                self._data = cast(ApiScalarType, {})
            elif self._metric_type == "LOG":
                self._data = cast(ApiScalarType, {})
            else:
                # 默认兜底到 scalar，实际上在实例化时被拦截
                self._data = cast(ApiScalarType, {})
        return cast(dict, self._data)

    @property
    def project_id(self) -> str:
        return self._project_id

    @property
    def run_id(self) -> str:
        return self._run_id

    @property
    def key(self) -> str:
        return self._key if self._key else ""

    @property
    def metric_type(self) -> str:
        return self._metric_type

    @property
    def metrics(self) -> List[Any]:
        return self._ensure_data().get("metrics", [])

    def _fetch_scalar(self) -> ApiScalarSeriesType:
        res = ApiScalarSeriesType(projectId=self.project_id, experimentId=self.run_id, key=self.key)
        # 1. 获取单指标列
        payload = {
            "projectId": self.project_id,
            "xType": "step",
            "range": [0, 0],
            "columns": [{"experimentId": self.run_id, "key": self.key}],
        }
        raw_resp = self._post("/house/metrics/scalar", data=payload)
        resp_list = (
            raw_resp.data if raw_resp.ok and isinstance(raw_resp.data, list) and len(raw_resp.data) > 0 else None
        )
        if resp_list is None:
            return res
        raw_data = resp_list[0]
        res["metrics"] = raw_data.get("metrics", {})
        # 2. 获取统计值列
        stat_resp = self._post("/house/metrics/scalar/value", data=payload)
        stat_list = (
            stat_resp.data if stat_resp.ok and isinstance(stat_resp.data, list) and len(stat_resp.data) > 0 else None
        )
        if stat_list is None:
            return res
        stat_data = stat_list[0]
        res["min"] = stat_data.get("min", {})
        res["max"] = stat_data.get("max", {})
        res["avg"] = stat_data.get("avg", {})
        res["median"] = stat_data.get("median", {})
        res["latest"] = stat_data.get("latest", {})
        return res

    def _fetch_media(self) -> Dict[str, Any]:
        return {}

    def _fetch_logs(self) -> Dict[str, Any]:
        return {}

    def export_csv(self) -> ApiResponseType:
        """
        导出列数据为 CSV。(同时支持 column 和 csv 导出)

        :return: ApiResponseType，成功时 data 包含临时下载 URL
        """
        resp = self._get(f"/experiment/{self._run_id}/column/csv", params={"key": self.key})
        if not resp.ok:
            return resp

        data = resp.data
        if isinstance(data, list) and data:
            url = data[0].get("url", "")
            return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))
        elif isinstance(data, dict):
            url = data.get("url", "")
            return ApiResponseType(ok=True, data=ApiColumnCsvExportType(url=url))

        return ApiResponseType(ok=False, errmsg="Invalid response format", data=None)

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
