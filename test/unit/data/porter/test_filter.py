"""
@author: cunyue
@file: test_filter.py
@time: 2025/7/19 21:26
@description: 测试过滤工具函数
"""

import pytest

from swanlab.data.porter.utils import filter_epoch, filter_metric, filter_column


@pytest.fixture
def mock_metrics():
    return {"key1": ("mock", "mock", None, 1)}


def test_returns_true_when_key_not_in_metrics(mock_metrics):
    assert filter_column("key2", mock_metrics) is True


def test_returns_false_when_key_in_metrics(mock_metrics):
    assert filter_column("key1", mock_metrics) is False


def test_returns_true_when_key_not_in_metrics_for_metric_filter(mock_metrics):
    assert filter_metric("key2", 2, mock_metrics) is True


def test_returns_true_when_step_greater_than_latest_step(mock_metrics):
    assert filter_metric("key1", 2, mock_metrics) is True


def test_returns_false_when_step_not_greater_than_latest_step(mock_metrics):
    assert filter_metric("key1", 1, mock_metrics) is False


def test_returns_true_when_epoch_greater_than_now_epoch():
    assert filter_epoch(5, 3) is True


def test_returns_false_when_epoch_not_greater_than_now_epoch():
    assert filter_epoch(3, 5) is False


def test_returns_true_when_now_epoch_is_none():
    assert filter_epoch(5, None) is True
