"""
@author: caddiesnew
@file: summary.py
@time: 2026/5/26
@description: Summary 实体类 — 实验标量指标概要统计查询
"""

from typing import Any, Dict, List, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.metric import ApiScalarSummaryItemType


class Summary(BaseEntity):
    """查询实验的标量指标概要统计数据（step/value + min/max/avg/median/stdDev）。"""

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        project_id: str,
        experiment_id: str,
        keys: Optional[List[str]] = None,
        root_pro_id: str = "",
        root_exp_id: str = "",
    ) -> None:
        super().__init__(ctx)
        self._project_id = project_id
        self._experiment_id = experiment_id
        self._keys = keys
        self._root_pro_id = root_pro_id
        self._root_exp_id = root_exp_id
        self._data: Optional[Dict[str, Any]] = None

    def _build_experiment_ref(self) -> Dict[str, str]:
        ref: Dict[str, str] = {
            "projectId": self._project_id,
            "experimentId": self._experiment_id,
        }
        if self._root_pro_id:
            ref["rootProId"] = self._root_pro_id
        if self._root_exp_id:
            ref["rootExpId"] = self._root_exp_id
        return ref

    def _ensure_data(self) -> Dict[str, ApiScalarSummaryItemType]:
        if self._data is None:
            payload: Dict[str, Any] = {"experiments": [self._build_experiment_ref()]}
            if self._keys is not None:
                payload["keys"] = self._keys
            resp = self._post("/house/metrics/summaries", data=payload)
            if resp.ok and isinstance(resp.data, dict):
                self._data = resp.data.get(self._experiment_id, {})
            else:
                self._data = {}
        assert self._data is not None
        return self._data

    @property
    def project_id(self) -> str:
        return self._project_id

    @property
    def experiment_id(self) -> str:
        return self._experiment_id

    def json(self) -> Dict[str, Any]:
        return dict(self._ensure_data())
