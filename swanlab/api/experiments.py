"""
@author: caddiesnew
@file: experiments.py
@time: 2026/4/20
@description: Experiments 迭代器 — 项目下的实验列表
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Union

from .base import BaseEntity
from .experiment import Experiment
from .utils import parse_column_type, to_camel_case

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


def _flatten_runs(runs: Union[List, Dict]) -> List:
    """展开分组后的实验数据，返回一个包含所有实验的列表。"""
    flat_runs = []
    items = runs.values() if isinstance(runs, dict) else [runs]
    for group in items:
        if isinstance(group, dict):
            flat_runs.extend(_flatten_runs(group))
        elif isinstance(group, list):
            flat_runs.extend(group)
    return flat_runs


class Experiments(BaseEntity):
    """
    项目下实验集合的迭代器。

    用法::

        for run in api.runs("username/project"):
            print(run.name)
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        path: str,
        filters: Optional[Dict[str, object]] = None,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._path = path
        self._filters = filters

    def __iter__(self) -> Iterator[Experiment]:
        parsed_filters = (
            [
                {
                    "key": to_camel_case(key) if parse_column_type(key) == "STABLE" else key.split(".", 1)[-1],
                    "active": True,
                    "value": [value],
                    "op": "EQ",
                    "type": parse_column_type(key),
                }
                for key, value in self._filters.items()
            ]
            if self._filters
            else []
        )
        resp = self._post(f"/project/{self._path}/runs/shows", data={"filters": parsed_filters})
        runs: List = []
        if isinstance(resp, list):
            runs = resp
        elif isinstance(resp, dict):
            runs = _flatten_runs(resp)

        for run_data in runs:
            yield Experiment(self._client, self._web_host, self._api_host, path=self._path, data=run_data)

    def to_dict(self) -> Dict[str, Any]:
        return {"path": self._path}
