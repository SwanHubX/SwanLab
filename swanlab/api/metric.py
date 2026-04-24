"""
@author: caddiesnew
@file: column.py
@time: 2026/4/20
@description: Column 实体类 — 实验列的查询与操作
"""

from typing import Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings import ApiColumnCsvExportType, ApiMetricTypeLiteral, ApiResponseType
from swanlab.api.typings.metric import ApiLogType, ApiMediaType, ApiMetricType, ApiScalarType
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

    def _ensure_data(self) -> ApiMetricType:
        if self._data is None:
            if self._metric_type == "SCALAR":
                self._data = cast(ApiScalarType, {})
            elif self._metric_type == "MEDIA":
                self._data = cast(ApiScalarType, {})
            elif self._metric_type == "LOG":
                self._data = cast(ApiScalarType, {})
            else:
                # 默认兜底到 scalar，实际上在实例化时被拦截
                self._data = cast(ApiScalarType, {})
        return self._data

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
        return []

    def _fetch_scalar(self):
        return

    def _fetch_media(self):
        pass

    def _fetch_logs(self):
        pass

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
