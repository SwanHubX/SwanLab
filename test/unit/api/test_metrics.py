"""
@author: Zhou QiYang
@file: test_metrics.py
@time: 2026/1/11 16:36
@description: 测试 Experiment.metrics() 方法，使用 MagicMock 和 monkeypatch 模拟网络请求
"""

from unittest.mock import patch, MagicMock

import pytest
import requests_mock

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
    分别模拟获取 csv 网址和文件内容
    """

    def __init__(self, metrics_data):
        self.metrics_data = metrics_data

    def __enter__(self):
        self.mock_get_metrics = patch('swanlab.api.experiment.thread.get_experiment_metrics').start()
        self.mock_get_metrics.side_effect = lambda client, expid, key: {'url': f'{get_host_api()}/{key}'}

        self.m = requests_mock.Mocker()
        self.m.start()
        for metric in self.metrics_data:
            self.m.get(f'{get_host_api()}/{metric[1]}', content=create_csv_data(*metric))
        return self

    def __exit__(self, *args):
        patch.stopall()
        self.m.stop()


def test_metrics_basic(experiment, metrics_data):
    """测试使用指定 keys 获取历史数据"""
    with MockSetup(metrics_data):
        result = experiment.metrics(keys=['loss', 'accuracy'])

    assert len(result) == 10
    assert 'loss' in result.columns
    assert 'accuracy' in result.columns


def test_metrics_with_x_axis(experiment, metrics_data):
    """测试使用 x_axis 参数"""
    with MockSetup(metrics_data):
        result = experiment.metrics(keys=['loss'], x_axis='accuracy')

    # x_axis 应该作为索引
    assert result.index.name == 'accuracy'


def test_metrics_with_sample(experiment, metrics_data):
    """测试使用 sample 参数限制返回行数"""
    with MockSetup(metrics_data):
        result = experiment.metrics(keys=['loss'], sample=5)

    # 只返回前 5 行
    assert len(result) == 5


def test_metrics_dict_mode(experiment, metrics_data):
    """测试 pandas=False 时返回 dict 格式"""
    with MockSetup(metrics_data):
        result = experiment.metrics(keys=['loss'], pandas=False)

    # 应该返回字典列表
    assert all(isinstance(item, dict) for item in result)


def test_full_metrics(experiment, metrics_data):
    """测试 keys 和 x_axis 都为 None 时调用 __full_metrics"""
    with MockSetup(metrics_data):
        result = experiment.metrics()

    assert len(result) == 10
    assert 'loss' in result.columns
    assert 'accuracy' in result.columns


@pytest.mark.parametrize("keys", ('invalid_keys', ['loss', 123, 'accuracy']))
def test_metrics_invalid_keys(experiment, metrics_data, keys):
    """测试 keys 参数类型错误的情况，返回空 DataFrame"""
    with MockSetup(metrics_data):
        result = experiment.metrics(keys=keys)
    assert len(result) == 0
