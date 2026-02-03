"""
@author: Zhou QiYang
@file: experiment.py
@time: 2026/1/11 16:36
@description: OpenApi 的单个实验对象
"""

from typing import List, Dict, Any

from swanlab.api.user import User
from swanlab.api.utils import Label, get_properties
from swanlab.core_python.api.type import RunType
from swanlab.core_python.client import Client
from swanlab.log import swanlog



class Experiment:
    def __init__(
        self, client: Client, *, data: RunType, path: str, web_host: str, login_user: str, line_count: int
    ) -> None:
        self._client = client
        self._data = data
        self._path = path
        self._web_host = web_host
        self._login_user = login_user
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
        return [Label(label['name']) for label in self._data['labels']]

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
        return User(client=self._client, login_user=self._login_user, username=self._data['user']['username'])

    @property
    def metric_keys(self) -> List[str]:
        """
        List of metric keys.
        """
        return list(self.summary.keys())

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

    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)

    def metrics(self, keys: List[str] = None, x_axis: str = None, sample: int = None, pandas: bool = True) -> Any:
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
        print(exp.metrics(keys=['loss'], sample=20, x_axis='t/accuracy'))

        Returns:
             t/accuracy    loss
        step
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
                "OpenApi requires pandas to implement the run.metrics(). Please install with 'pip install pandas'."
            )
        
        if keys is None:
            raise ValueError("keys cannot be None")
        
        x_axis_state = x_axis is not None and x_axis != "step"

        if isinstance(keys, str):
            keys = [keys]
        if x_axis_state:
            keys += [x_axis]

        # 去重 keys
        keys = list(set(keys))
        dfs = []
        prefix = ""
        for idx, key in enumerate(keys):
            resp = self._client.get(f"/experiment/{self.id}/column/csv", params={"key": key})

            url:str = resp[0].get("url", "")
            df = pd.read_csv(url, index_col=0)

            if idx == 0:
                # 从第一列名提取 prefix，例如 "t0707-02:17-loss_step" 中提取 "t0707-02:17-"
                first_col = df.columns[0]
                suffix = f"{key}_"
                if suffix in first_col:
                    prefix = first_col.split(suffix)[0]  # 结果为 "t0707-02:17-"
                else:
                    prefix = ""

            if prefix:
                df.columns = [
                    col[len(prefix):].removesuffix("_step") if col.startswith(prefix) else col.removesuffix("_step")
                    for col in df.columns
                ]
            else:
                df.columns = [col.removesuffix("_step") for col in df.columns]
            
            dfs.append(df)

        # 拼接整张表
        result_df = dfs[0]
        if len(dfs) > 1:
            for df in dfs[1:]:
                result_df = result_df.join(df, how='outer').sort_index()

        # 如果有 x_axis，进行特殊处理
        if x_axis_state:
            # 去掉所有带 _timestamp 后缀的列
            timestamp_cols = [col for col in result_df.columns if col.endswith("_timestamp")]
            result_df = result_df.drop(columns=timestamp_cols)
            
            # 确保 x_axis 列存在
            if x_axis not in result_df.columns:
                raise ValueError(f"x_axis '{x_axis}' not found in the result DataFrame")
            
            # 将 x_axis 列放到第一列
            cols = [x_axis] + [col for col in result_df.columns if col != x_axis]
            result_df = result_df[cols]
            result_df = result_df[result_df[x_axis].notna()]
        
        if sample is not None:
            result_df = result_df.head(sample)

        return result_df


__all__ = ['Experiment']
