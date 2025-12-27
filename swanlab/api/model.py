"""
@author: Zhou Qiyang
@file: model.py
@time: 2025/12/18 20:10
@description: OpenApi查询结果将以对象返回，并且对后端的返回字段进行一些筛选
"""

from typing import List, Dict, Optional

from swanlab.core_python import Client
from swanlab.core_python.api.experiment import get_project_experiments, get_experiment_metrics
from swanlab.core_python.api.project import get_workspace_projects
from swanlab.core_python.api.type import ProjectType, ProjectLabelType, ProjResponseType, UserType, RunType
from .utils import flatten_runs


class Label:
    """
    Project label object
    you can get the label name by str(label)
    """

    def __init__(self, data: ProjectLabelType):
        self._data = data

    @property
    def name(self):
        """
        Label name.
        """
        return self._data['name']

    def __str__(self):
        return str(self.name)


class User:
    def __init__(self, data: UserType):
        self._data = data

    @property
    def name(self) -> str:
        return self._data['name']

    @property
    def username(self) -> str:
        return self._data['username']


class Project:
    """
    Representing a single project with some of its properties.
    """

    def __init__(self, data: ProjectType, web_host: str):
        self._data = data
        self._web_host = web_host

    @property
    def name(self) -> str:
        """
        Project name.
        """
        return self._data['name']

    @property
    def path(self) -> str:
        """
        Project path in the format 'username/project-name'.
        """
        return self._data['path']

    @property
    def url(self) -> str:
        """
        Full URL to access the project.
        """
        return f"{self._web_host}/@{self._data['path']}"

    @property
    def description(self) -> str:
        """
        Project description.
        """
        return self._data['description']

    @property
    def visibility(self) -> str:
        """
        Project visibility, either 'PUBLIC' or 'PRIVATE'.
        """
        return self._data['visibility']

    @property
    def created_at(self) -> str:
        """
        Project creation timestamp
        """
        return self._data['createdAt']

    @property
    def updated_at(self) -> str:
        """
        Project last update timestamp
        """
        return self._data['updatedAt']

    @property
    def workspace(self) -> str:
        """
        Project workspace name.
        """
        return self._data["group"]["username"]

    @property
    def labels(self) -> List[Label]:
        """
        List of Label attached to this project.
        """
        return [Label(label) for label in self._data['projectLabels']]

    @property
    def count(self) -> Dict[str, int]:
        """
        Project statistics dictionary containing:
        experiments, contributors, children, collaborators, runningExps.
        """
        return self._data['_count']

    @property
    def __dict__(self) -> Dict[str, object]:
        """
        Return a dictionary containing all @property fields.
        """
        result = {}
        cls = type(self)
        for attr_name in dir(cls):
            if attr_name.startswith('_'):
                continue
            attr = getattr(cls, attr_name, None)
            if isinstance(attr, property):
                result[attr_name] = self.__getattribute__(attr_name)
        return result


class Experiment:
    def __init__(self, data: RunType, client: Client, path: str, web_host: str, line_count: int):
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

    @property
    def __dict__(self) -> Dict[str, object]:
        """
        Return a dictionary containing all @property fields.
        """
        result = {}
        cls = type(self)
        for attr_name in dir(cls):
            if attr_name.startswith('_'):
                continue
            attr = getattr(cls, attr_name, None)
            if isinstance(attr, property):
                result[attr_name] = self.__getattribute__(attr_name)
        return result

    def history(self, keys: List[str], x_axis: str = None, sample: int = None, pandas: bool = True):
        """
        Get specific metric data of the experiment.
        :param keys: List of metric keys to obtain.
        :param x_axis: The metric to be used as x-axis. If None, 'step' will be used as the x-axis.
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
        df = pd.DataFrame()

        # x轴不为空时，将step替换为x_axis的指标数据
        if x_axis is not None:
            csv_df = get_experiment_metrics(self._client, expid=self.id, key=x_axis)
            df[x_axis] = csv_df.iloc[:, 1]
        for key in keys:
            csv_df = get_experiment_metrics(self._client, expid=self.id, key=key)
            if csv_df is None:
                continue
            if key == keys[0] and x_axis is None:
                df['step'] = df[key] = csv_df.iloc[:, 0]
            df[key] = csv_df.iloc[:, 1]

        # 截取前sample行
        if sample is not None:
            df = df.head(sample)
        return df if pandas else df.to_dict(orient='records')

    def scan_history(self, min_step: int = 0, max_step: int = 1000):
        """
        Get all metric keys' data of the experiment in a range of steps.
        :param min_step: The minimum step to get.
        :param max_step: The maximum step to get.
        :return: Metric data in dict.
        """
        return self.history(keys=self.metric_keys, pandas=False)[min_step:max_step]


class Projects:
    """
    Container for a collection of Project objects.
    You can iterate over the projects by for-in loop.
    """

    def __init__(
        self,
        client: Client,
        web_host: str,
        workspace: str,
        sort: Optional[List[str]] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ):
        self._client = client
        self._web_host = web_host
        self._workspace = workspace
        self._sort = sort
        self._search = search
        self._detail = detail

    def __iter__(self):
        # 按用户遍历情况获取项目信息
        cur_page = 0
        page_size = 20
        while True:
            cur_page += 1
            projects_info: ProjResponseType = get_workspace_projects(
                self._client,
                workspace=self._workspace,
                page=cur_page,
                size=page_size,
                sort=self._sort,
                search=self._search,
                detail=self._detail,
            )
            if cur_page * page_size >= projects_info['total']:
                break

        yield from iter(Project(project, self._web_host) for project in projects_info['list'])


class Experiments:
    """
    Container for a collection of Experiment objects.
    You can iterate over the experiments by for-in loop.
    """

    def __init__(self, client: Client, path: str, web_host: str, filters: Dict[str, object] = None):
        if len(path.split('/')) != 2:
            raise ValueError(f"User's {path} is invaded. Correct path should be like 'username/project'")
        self._client = client
        self._path = path
        self._web_host = web_host
        self._filters = filters

    def __iter__(self):
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
