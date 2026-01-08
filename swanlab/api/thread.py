"""
@author: Zhou QiYang
@file: thread.py
@time: 2025/12/30 15:08
@description: 用于api并发请求的封装类
"""

from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from typing import List, Any, TYPE_CHECKING

import requests

if TYPE_CHECKING:
    from swanlab.core_python.client import Client

from swanlab.core_python.api.experiment import get_experiment_metrics
from swanlab.log import swanlog


class HistoryPool:

    def __init__(self, client: "Client", expid: str, *, keys: List[str], x_axis: str = None, num_threads: int = 10):
        try:
            import pandas as pd
        except ImportError:
            raise TypeError(
                "OpenApi requires pandas to init the HistoryPool. Please install with 'pip install pandas'."
            )

        self._client = client
        self._expid = expid
        self._keys = keys
        self._x_axis = x_axis
        if self._x_axis is not None:
            self._keys = [self._x_axis] + [k for k in self._keys if k != self._x_axis]
        self._num_threads = num_threads

        # 使用 _results 字典收集每个 key 的 DataFrame，最后统一按顺序合并到 _history
        self._executor = ThreadPoolExecutor(max_workers=self._num_threads)
        self._futures = []
        self._results = dict()
        self._history = pd.DataFrame()

    def _task(self, key: str):
        """
        处理单个key，获取对应csv
        """
        import pandas as pd

        try:
            csv_df = pd.DataFrame()
            resp = get_experiment_metrics(self._client, expid=self._expid, key=key)
            # 从返回网址中解析csv内容
            with requests.get(resp['url']) as response:
                csv_df = pd.read_csv(BytesIO(response.content))
            return key, csv_df
        except Exception as e:
            swanlog.warning(f'Error processing key {key} in experiment {self._expid}: {e}')
            return key, pd.DataFrame()

    def execute(self) -> Any:
        if not self._keys:
            return self._history

        # 将所有key提交到线程池
        for key in self._keys:
            future = self._executor.submit(self._task, key)
            self._futures.append((key, future))

        # 等待所有任务完成并收集结果
        for key, future in self._futures:
            try:
                result_key, csv_df = future.result()
                self._results[result_key] = csv_df
            except Exception as e:
                swanlog.warning(f'Error getting result for key {key} in experiment {self._expid}: {e}')
        self._executor.shutdown(wait=True)

        # 按照 keys 的顺序统一合并
        for key in self._keys:
            if key not in self._results:
                continue
            key_df = self._results[key]
            step_col, value_col = key_df.columns[:2]  # step 列, 指标值列

            # 将 step 设为索引，其后基于索引自动对齐
            if self._history.empty:
                self._history = key_df.set_index(step_col)
            else:
                self._history[value_col] = key_df.set_index(step_col)[value_col]

        # 若指定x轴，重置索引
        if self._x_axis is not None:
            self._history = self._history.reset_index().iloc[:, 1:]
            self._history = self._history.set_index(self._history.columns[0])
        else:
            self._history.rename(columns={'step': '_step'}, inplace=True)
        return self._history
