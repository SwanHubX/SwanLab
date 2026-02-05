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


class Profile:
    """Experiment profile containing config, metadata, requirements, and conda info."""

    def __init__(self, data: Dict):
        self._data = data

    @property
    def config(self) -> Dict:
        """Experiment configuration."""
        return self._data.get('config', {})

    @property
    def metadata(self) -> Dict:
        """Experiment metadata."""
        return self._data.get('metadata', {})

    @property
    def requirements(self) -> str:
        """Python requirements."""
        return self._data.get('requirements', '')

    @property
    def conda(self) -> str:
        """Conda environment."""
        return self._data.get('conda', '')

    @property
    def scalar(self) -> Dict:
        """Experiment scalar metrics."""
        return self._data.get('scalar', {})


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
        return self._data.get('name', '')

    @property
    def id(self) -> str:
        """
        Experiment CUID.
        """
        return self._data.get('cuid', '')

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
        return self._data.get('createdAt', '')

    @property
    def finished_at(self) -> str:
        """
        Experiment finished timestamp
        """
        return self._data.get('finishedAt', '')

    @property
    def profile(self) -> Profile:
        """
        Experiment profile containing config, metadata, requirements, and conda.
        """
        return Profile(self._data.get('profile', {}))

    @property
    def show(self) -> bool:
        """
        Whether the experiment is visible.
        """
        return self._data.get('show', True)

    @property
    def description(self) -> str:
        """
        Experiment description.
        """
        return self._data.get('description', '')

    @property
    def labels(self) -> List[Label]:
        """
        List of Label attached to this experiment.
        """
        return [Label(label['name']) for label in self._data.get('labels', [])]

    @property
    def state(self) -> str:
        """
        Experiment state.
        """
        return self._data.get('state', '')

    @property
    def group(self) -> str:
        """
        Experiment group.
        """
        return self._data.get('cluster', '')

    @property
    def job(self) -> str:
        """
        Experiment job type.
        """
        return self._data.get('job', '')

    @property
    def user(self) -> User:
        """
        Experiment user.
        """
        username = self._data.get('user', {}).get('username', '')
        return User(client=self._client, login_user=self._login_user, username=username)

    @property
    def metric_keys(self) -> List[str]:
        """
        List of metric keys.
        """
        return list(self.profile.scalar.keys())

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
        return self._data.get('rootExpId', '')

    @property
    def root_pro_id(self) -> str:
        """
        Root project cuid. If the experiment is a root experiment, it will be None.
        """
        return self._data.get('rootProId', '')

    def metrics(self, keys: List[str] = None, x_axis: str = None, sample: int = None, pandas: bool = True) -> Any:
        """
        Get metric data from the experiment.

        Args:
            keys: List of metric keys to fetch. Required. A single string is also accepted.
            x_axis: Metric to use as x-axis. Defaults to 'step'.
            sample: Number of rows to return from the start.
            pandas: Reserved parameter (always returns DataFrame).

        Returns:
            pandas.DataFrame with metric data.

        Example:
            >>> exp.metrics(keys=['loss', 'accuracy'], sample=20, x_axis='t/accuracy')
                 t/accuracy    loss
            step
            0    0.310770   0.525776
            1    0.642817   0.479186
            ...
        """
        try:
            import pandas as pd
        except ImportError:
            raise TypeError("pandas is required for metrics(). Install with 'pip install pandas'.")

        # Normalize keys: must be a non-empty list of strings
        if keys is None:
            swanlog.warning('keys cannot be None')
            return pd.DataFrame()
        if not isinstance(keys, list):
            swanlog.warning('keys must be a list')
            return pd.DataFrame()
        if not keys:
            swanlog.warning('keys cannot be empty')
            return pd.DataFrame()
        if not all(isinstance(k, str) for k in keys):
            swanlog.warning('keys must be a list of strings')
            return pd.DataFrame()

        # Determine if x_axis needs to be included
        use_x_axis = x_axis is not None and x_axis != "step"
        if use_x_axis:
            keys.append(x_axis)

        # Fetch and process each metric CSV
        dfs = []
        prefix = ""
        for idx, key in enumerate(keys):
            resp = self._client.get(f"/experiment/{self.id}/column/csv", params={"key": key})
            url = resp[0].get("url", "")
            df = pd.read_csv(url, index_col=0)

            # Extract prefix from first column (e.g., "t0707-02:17-loss_step" → "t0707-02:17-")
            if idx == 0:
                first_col = df.columns[0]
                suffix = f"{key}_"
                prefix = first_col.split(suffix)[0] if suffix in first_col else ""

            # Strip "_step" suffix from column names (Python 3.8 compatible)
            def strip_suffix(col, suffix="_step"):
                return col[:-len(suffix)] if col.endswith(suffix) else col

            # Apply prefix removal and suffix stripping
            df.columns = [
                strip_suffix(col[len(prefix):]) if prefix and col.startswith(prefix) else strip_suffix(col)
                for col in df.columns
            ]
            dfs.append(df)

        # Merge all DataFrames
        result_df = dfs[0].join(dfs[1:], how='outer') if len(dfs) > 1 else dfs[0]
        result_df = result_df.sort_index()

        # Handle x_axis: drop timestamp columns, reorder, filter nulls
        if use_x_axis:
            result_df = result_df.drop(columns=[c for c in result_df.columns if c.endswith("_timestamp")], errors='ignore')
            if x_axis not in result_df.columns:
                raise ValueError(f"x_axis '{x_axis}' not found in result DataFrame")
            cols = [x_axis] + [c for c in result_df.columns if c != x_axis]
            result_df = result_df[cols].dropna(subset=[x_axis])

        # Apply sample limit
        if sample is not None:
            result_df = result_df.head(sample)

        return result_df
    
    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)


__all__ = ['Experiment', 'Profile']
