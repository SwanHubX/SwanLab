"""
@author: Zhou QiYang
@file: thread.py
@time: 2025/12/30 15:08
@description: OpenApi 的实验对象迭代器
"""

from typing import List, Dict, Iterator

from swanlab.api.experiment import Experiment
from swanlab.api.utils import ApiBase
from swanlab.core_python.api.experiment import get_project_experiments
from swanlab.core_python.api.type import RunType
from swanlab.core_python.auth.providers.api_key import LoginInfo
from swanlab.core_python.client import Client


def flatten_runs(runs: Dict) -> List:
    """
    展开分组后的实验数据，返回一个包含所有实验的列表
    """
    flat_runs = []
    for group in runs.values():
        if isinstance(group, Dict):
            flat_runs.extend(flatten_runs(group))
        else:
            flat_runs.extend(group)
    return flat_runs


class Experiments(ApiBase):
    """
    Container for a collection of Experiment objects.
    You can iterate over the experiments by for-in loop.
    """

    def __init__(self, client: Client, *, path: str, login_info: LoginInfo, filters: Dict[str, object] = None) -> None:
        if len(path.split('/')) != 2:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project'")
        super().__init__(client)
        self._path = path
        self._web_host = login_info.web_host
        self._login_user = login_info.username
        self._filters = filters

    def __iter__(self) -> Iterator[Experiment]:
        # TODO: 完善filter的功能（正则、条件判断）
        resp = get_project_experiments(self._client, path=self._path, filters=self._filters)
        runs: List[RunType] = []
        if isinstance(resp, List):
            runs = resp
        # 分组时需展平实验数据
        elif isinstance(resp, Dict):
            runs = flatten_runs(resp)
        line_count = len(runs)
        yield from iter(
            Experiment(
                self._client,
                data=run,
                path=self._path,
                web_host=self._web_host,
                login_user=self._login_user,
                line_count=line_count,
            )
            for run in runs
        )


__all__ = ["Experiments"]
