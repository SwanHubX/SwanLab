"""
@author: Zhou QiYang
@file: test_metrics.py
@time: 2026/1/11 16:36
@description: 测试 Experiment.metrics() 方法，使用 MagicMock 和 monkeypatch 模拟网络请求
"""

from unittest.mock import patch, MagicMock

import pytest
import pandas as pd

from swanlab.api.experiment import Experiment
from swanlab.core_python.client import Client
from swanlab.package import get_host_web, get_host_api
from utils import create_csv_data, create_run_type_data


@pytest.fixture
def experiment():
    """创建 Experiment 实例"""
    data = create_run_type_data()
    return Experiment(
        MagicMock(spec=Client),
        data=data,
        path='test_user/test_project/test_exp',
        web_host=get_host_web(),
        login_user='test_user',
        line_count=100,
    )


@pytest.fixture
def metrics_data():
    return [
        (list(range(10)), 'loss', [0.5 - i * 0.05 for i in range(10)]),
        (list(range(10)), 'accuracy', [0.5 + i * 0.05 for i in range(10)]),
    ]


class MockSetup:
    """
    模拟网络请求
    直接 mock client.get 返回 URL，然后 mock pd.read_csv 返回 DataFrame
    """

    def __init__(self, metrics_data, experiment):
        self.metrics_data = metrics_data
        self.experiment = experiment
        # Create a lookup dict for metrics data
        self._metric_lookup = {m[1]: m for m in metrics_data}

    def _create_df(self, key):
        """Helper to create DataFrame from metric data"""
        if key not in self._metric_lookup:
            return pd.DataFrame()
        step_values, metric_name, metric_values = self._metric_lookup[key]
        return pd.DataFrame({
            'step': step_values,
            f'{metric_name}_step': metric_values
        }).set_index('step')

    def __enter__(self):
        # Mock client.get to return CSV URL
        self.mock_get = patch.object(self.experiment._client, 'get').start()
        self.mock_get.side_effect = lambda path, params: [{'url': f'{get_host_api()}/{params["key"]}'}]

        # Mock pd.read_csv to return DataFrame directly
        self.mock_read_csv = patch('pandas.read_csv').start()
        self.mock_read_csv.side_effect = lambda url, index_col: self._create_df(url.split('/')[-1])
        return self

    def __exit__(self, *args):
        patch.stopall()


def test_metrics_basic(experiment, metrics_data):
    """测试使用指定 keys 获取历史数据"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics(keys=['loss', 'accuracy'])

    assert len(result) == 10
    assert 'loss' in result.columns
    assert 'accuracy' in result.columns


def test_metrics_with_x_axis(experiment, metrics_data):
    """测试使用 x_axis 参数"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics(keys=['loss'], x_axis='accuracy')

    # x_axis 列应该在第一列
    assert result.columns[0] == 'accuracy'


def test_metrics_with_sample(experiment, metrics_data):
    """测试使用 sample 参数限制返回行数"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics(keys=['loss'], sample=5)

    # 只返回前 5 行
    assert len(result) == 5


def test_metrics_dict_mode(experiment, metrics_data):
    """测试 pandas=False 时返回 DataFrame（当前实现只支持 DataFrame）"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics(keys=['loss'], pandas=False)

    # 当前实现始终返回 DataFrame
    assert isinstance(result, pd.DataFrame)


def test_full_metrics(experiment, metrics_data):
    """测试 keys=None 时返回空 DataFrame（当前实现不支持）"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics()

    # keys=None 时返回空 DataFrame
    assert len(result) == 0


@pytest.mark.parametrize("keys", ('invalid_keys', ['loss', 123, 'accuracy']))
def test_metrics_invalid_keys(experiment, metrics_data, keys):
    """测试 keys 参数类型错误的情况，返回空 DataFrame"""
    with MockSetup(metrics_data, experiment):
        result = experiment.metrics(keys=keys)
    assert len(result) == 0
