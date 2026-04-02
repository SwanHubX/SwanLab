"""
@author: cunyue
@file: test_upload.py
@time: 2025/11/9 18:55
@description: 测试上传
"""

import time
from typing import List
from unittest.mock import MagicMock, patch

import pytest

from swanlab.core_python.uploader.batch import MetricDict, trace_metrics
from swanlab.core_python.uploader.model import ScalarModel
from swanlab.core_python.uploader.thread.log_collector import LogCollectorTask
from swanlab.core_python.uploader.thread.task_types import UploadType
from swanlab.core_python.uploader.upload import dedupe_metrics_by_key_step


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


def test_dedupe_metrics_by_key_step_keeps_highest_epoch_duplicate():
    metrics = [
        ScalarModel({"data": 3.0}, "loss", 51, 52),
        ScalarModel({"data": 3.0}, "loss", 51, 52),
        ScalarModel({"data": 2.0}, "loss", 50, 53),
        ScalarModel({"data": 1.0}, "loss", 50, 51),
    ]

    deduped = dedupe_metrics_by_key_step(metrics)

    assert len(deduped) == 2
    assert {
        (metric.key, metric.step, metric.epoch, metric.metric["data"])
        for metric in deduped
    } == {
        ("loss", 50, 53, 2.0),
        ("loss", 51, 52, 3.0),
    }


def test_log_collector_upload_dedupes_duplicate_scalar_steps():
    first = ScalarModel({"data": 1.0}, "loss", 50, 51)
    second = ScalarModel({"data": 2.0}, "loss", 50, 51)
    captured = {}

    def fake_upload(metrics, upload_callback=None):
        captured["metrics"] = metrics
        return None, None

    collector = LogCollectorTask()
    collector.container = [(UploadType.SCALAR_METRIC, [first]), (UploadType.SCALAR_METRIC, [second])]

    with patch.dict(UploadType.SCALAR_METRIC.value, {"upload": fake_upload}, clear=False):
        collector.upload()

    assert len(captured["metrics"]) == 1
    assert captured["metrics"][0] is second
    assert all(task_type != UploadType.SCALAR_METRIC for task_type, _ in collector.container)
