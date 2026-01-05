"""
@author: Zhou QiYang
@file: experiment.py
@time: 2026/1/5 17:58
@description: OpenApi 中的实验对象
"""

from typing import TYPE_CHECKING, List, Dict, Any, Iterator

if TYPE_CHECKING:
    from swanlab.core_python.client import Client

from swanlab.log import swanlog
from swanlab.core_python.api.experiment import get_project_experiments
from swanlab.core_python.api.type import RunType
from swanlab.api.thread import HistoryPool
from swanlab.api.utils import flatten_runs

from .base import ApiBase, Label, User


class Experiment(ApiBase):
    def __init__(self, data: RunType, client: "Client", path: str, web_host: str, line_count: int) -> None:
        self._data = data
        self._client = client
        self._path = path
        self._web_host = web_host
        self._line_count = line_count

    @property
    def name(self) -> str:
        """
        Experiment name.
        """
        return self._data['name']

    @property
    def id(self) -> str:
        """
        Experiment CUID.
        """
        return self._data['cuid']

    @property
    def url(self) -> str:
        """
        Full URL to access the experiment.
        """
        return f"{self._web_host}/@{self._path}/runs/{self.id}/chart"

    @property
    def created_at(self) -> str:
        """
        Experiment creation timestamp
        """
        return self._data['createdAt']

    @property
    def description(self) -> str:
        """
        Experiment description.
        """
        return self._data['description']

    @property
    def labels(self) -> List[Label]:
        """
        List of Label attached to this experiment.
        """
        return [Label(label) for label in self._data['labels']]

    @property
    def config(self) -> Dict[str, object]:
        """
        Experiment configuration. Can be used as filter in the format of 'config.<key>'
        """
        return self._data['profile']['config']

    @property
    def summary(self) -> Dict[str, object]:
        """
        Experiment metrics data. Can be used as filter in the format of 'summary.<key>'
        """
        return self._data['profile']['scalar']

    @property
    def state(self) -> str:
        """
        Experiment state.
        """
        return self._data['state']

    @property
    def group(self) -> str:
        """
        Experiment group.
        """
        return self._data['cluster']

    @property
    def job(self) -> str:
        """
        Experiment job type.
        """
        return self._data['job']

    @property
    def user(self) -> User:
        """
        Experiment user.
        """
        return User(self._data['user'])

    @property
    def metric_keys(self) -> List[str]:
        """
        List of metric keys.
        """
        summary_keys = self.summary.keys()
        return list(summary_keys)

    @property
    def history_line_count(self) -> int:
        """
        The number of historical experiments in this project.
        """
        return self._line_count

    @property
    def root_exp_id(self) -> str:
        """
        Root experiment cuid. If the experiment is a root experiment, it will be None.
        """
        return self._data['rootExpId']

    @property
    def root_pro_id(self) -> str:
        """
        Root project cuid. If the experiment is a root experiment, it will be None.
        """
        return self._data['rootProId']

    def __full_history(self) -> Any:
        """
        Get all metric keys' data of the experiment with timestamp.
        """
        try:
            import pandas as pd
        except ImportError:
            raise TypeError(
                "OpenApi requires pandas to implement the run.history(). Please install with 'pip install pandas'."
            )

        df = pd.DataFrame()
        if len(self.metric_keys) >= 1:
            pool = HistoryPool(self._client, self.id, self.metric_keys)
            pool.start()
            df = pool.wait_completion()

        return df

    def history(self, keys: List[str] = None, x_axis: str = None, sample: int = None, pandas: bool = True) -> Any:
        """
        Get specific metric data of the experiment.
        :param keys: List of metric keys to obtain. If None, all metrics keys will be used.
        :param x_axis: The metric to be used as x-axis. If None, '_step' will be used as the x-axis.
        :param sample: Number of rows to select from the beginning.
        :param pandas: Whether to return a pandas DataFrame. If False, returns dict format: {key: [values], ...}
        :return: Metric data.

        Example:
        ```python
        api = swanlab.OpenApi()
        exp = api.run(path="username/project/expid")  # You can get expid from api.runs()
        print(exp.history(keys=['loss'], sample=20, x_axis='t/accuracy'))

        Returns:
            t/accuracy    loss
        0    0.310770   0.525776
        1    0.642817   0.479186
        2    0.646031   0.362428
        3    0.608820   0.230555
        ...
        19   0.791999   0.180106
        ```
        """
        try:
            import pandas as pd
        except ImportError:
            raise TypeError(
                "OpenApi requires pandas to implement the run.history(). Please install with 'pip install pandas'."
            )

        if keys is not None and not isinstance(keys, list):
            swanlog.warning('keys must be specified as a list')
            return pd.DataFrame()
        elif keys is not None and len(keys) and not all(isinstance(k, str) for k in keys):
            swanlog.warning('keys must be a list of string')
            return pd.DataFrame()

        # 使用 merge 按 step 对齐不同指标的数据
        df = pd.DataFrame()

        # 使用线程池并发获取所有的key的指标数据
        if keys is not None:
            keys = [x_axis] + [k for k in keys if k != x_axis]
            pool = HistoryPool(self._client, self.id, keys)
            pool.start()
            df = pool.wait_completion()

        # x轴不为空时，将step替换为x_axis的指标数据
        if x_axis is not None:
            df = df.drop(columns=['_step'])
        # x轴与keys都未指定时，返回带时间戳的所有指标数据
        elif keys is None:
            df = self.__full_history()

        # 按 x 轴排序
        if x_axis in df.columns:
            df = df.sort_values(by=x_axis).reset_index(drop=True)

        # 截取前sample行
        if sample is not None:
            df = df.head(sample)

        return df if pandas else df.to_dict(orient='records')


class Experiments(ApiBase):
    """
    Container for a collection of Experiment objects.
    You can iterate over the experiments by for-in loop.
    """

    def __init__(self, client: "Client", path: str, web_host: str, filters: Dict[str, object] = None) -> None:
        if len(path.split('/')) != 2:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project'")
        self._client = client
        self._path = path
        self._web_host = web_host
        self._filters = filters

    def __iter__(self) -> Iterator[Experiment]:
        # todo: 完善filter的功能（正则、条件判断）
        resp = get_project_experiments(self._client, path=self._path, filters=self._filters)
        runs: List[RunType] = []
        if isinstance(resp, List):
            runs = resp
        # 分组时需展平实验数据
        elif isinstance(resp, Dict):
            runs = flatten_runs(resp)
        line_count = len(runs)
        yield from iter(Experiment(run, self._client, self._path, self._web_host, line_count) for run in runs)
