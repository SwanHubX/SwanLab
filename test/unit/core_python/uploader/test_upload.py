"""
@author: cunyue
@file: test_upload.py
@time: 2025/11/9 18:55
@description: 测试上传
"""

import time
from typing import List
from unittest.mock import patch, MagicMock

import pytest

from swanlab.core_python.uploader.batch import trace_metrics, MetricDict


@pytest.fixture
def mock_client():
    client = MagicMock()
    client.pending = False
    return client


def mock_metrics(metrics: List[dict]) -> MetricDict:
    return {
        "projectId": "proj_123",
        "experimentId": "exp_456",
        "type": "scalar",
        "metrics": metrics,
        "flagId": None,
    }


def test_trace_metrics_uploads_when_client_is_pending(mock_client):
    with patch("swanlab.core_python.uploader.batch.get_client", return_value=mock_client):
        mock_client.pending = True
        mock_client.post.return_value = (None, MagicMock(status_code=200))
        trace_metrics("/test/url", data=mock_metrics([{"key": "value"}]))
        # 即使client处于pending状态，也会尝试上传一次
        mock_client.post.assert_called_once()


def test_trace_metrics_uploads_all_metrics_in_batches(mock_client):
    with patch("swanlab.core_python.uploader.batch.get_client", return_value=mock_client):
        mock_client.post.return_value = (None, MagicMock(status_code=200))
        start = time.time()
        trace_metrics("/test/url", data=mock_metrics([{"key": "value"}] * 2500), per_request_len=1000)
        assert mock_client.post.call_count == 3
        end = time.time()
        assert end - start > 3
