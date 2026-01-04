"""
@author: Zhou QiYang
@file: thread.py
@time: 2025/12/30 15:08
@description: 用于api并发请求的封装类
"""

import queue
import threading
from io import BytesIO
from typing import List, Any, TYPE_CHECKING

import requests

if TYPE_CHECKING:
    from swanlab.core_python.client import Client

from swanlab.core_python.api.experiment import get_experiment_metrics
from swanlab.error import ApiError
from swanlab.log import swanlog


class HistoryPool:

    def __init__(self, client: "Client", expid: str, keys: List[str], num_threads: int = 10):
        try:
            import pandas as pd
        except ImportError:
            raise TypeError(
                "OpenApi requires pandas to init the HistoryPool. Please install with 'pip install pandas'."
            )

        self._client = client
        self._expid = expid
        self._keys = keys
        self._num_threads = min(num_threads, len(keys)) if keys else num_threads
        self._task_queue = queue.Queue()
        self._threads = []
        # 使用 _results 字典收集每个 key 的 DataFrame，最后统一 merge 到 _history
        self._history = pd.DataFrame()
        self._results = dict()
        self._results_lock = threading.Lock()

    def start(self):
        """
        启动工作线程
        """
        if not self._keys:
            return

        # 将所有key放入队列
        for key in self._keys:
            self._task_queue.put(key)

        # 创建固定数量的工作线程
        for _ in range(self._num_threads):
            thread = threading.Thread(target=self._worker)
            thread.daemon = True
            self._threads.append(thread)
            thread.start()

    def _worker(self):
        """
        工作线程函数，从队列中取 key ，获取对应csv添加到 _result 中
        """
        import pandas as pd

        while True:
            key = self._task_queue.get()
            if key is None:  # 结束信号
                self._task_queue.task_done()
                break

            try:
                csv_df = pd.DataFrame()
                resp = get_experiment_metrics(self._client, expid=self._expid, key=key)
                try:
                    # 从返回网址中解析csv内容
                    with requests.get(resp['url']) as response:
                        csv_df = pd.read_csv(BytesIO(response.content))
                        csv_df = csv_df.drop(csv_df.columns[-1], axis=1)
                except ApiError:
                    swanlog.warning(f'key {key} does not exist in experiment: {self._expid}')
                    pass
                # 只存储结果，不进行 merge
                with self._results_lock:
                    self._results[key] = csv_df
            except Exception as e:
                swanlog.warning(f'Error processing key {key} in experiment {self._expid}: {e}')
            finally:
                self._task_queue.task_done()

    def wait_completion(self) -> Any:
        import pandas as pd

        # 等待所有任务完成，发送结束信号给所有工作线程
        self._task_queue.join()
        for _ in range(self._num_threads):
            self._task_queue.put(None)
        for thread in self._threads:
            thread.join()

        # 按照 keys 的顺序统一合并
        for key in self._keys:
            if key not in self._results:
                continue
            key_df = self._results[key]
            if self._history.empty:
                self._history = key_df
            else:
                self._history = pd.merge(self._history, key_df, on='step', how='outer')

        self._history.rename(columns={'step': '_step'}, inplace=True)
        return self._history
