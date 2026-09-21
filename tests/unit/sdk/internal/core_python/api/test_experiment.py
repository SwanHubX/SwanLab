from datetime import datetime, timezone
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.run.v1.run_pb2 import RUN_STATE_ABORTED, RUN_STATE_FINISHED
from swanlab.sdk.internal.core_python.api import experiment as experiment_api


def _created_at() -> Timestamp:
    ts = Timestamp()
    ts.FromDatetime(datetime(2024, 8, 1))
    return ts


def _post_resp(data: dict, status_code: int) -> SimpleNamespace:
    return SimpleNamespace(data=data, raw=SimpleNamespace(status_code=status_code))


def _experiment_data(**extra) -> dict:
    data = {"cuid": "experiment-id", "slug": "run-id", "name": "test-run"}
    data.update(extra)
    return data


def _call_create_or_resume():
    return experiment_api.create_or_resume_experiment(
        "alice",
        "demo",
        name="test-run",
        resume="allow",
        color="#abcdef",
        description=None,
        job_type=None,
        group=None,
        tags=None,
        created_at=_created_at(),
    )


def test_create_new_experiment_reports_current_created_at(monkeypatch):
    # createdAt 上报的是指标开始上报的时刻（当前时间），而非本地记录的创建时间
    post = MagicMock(return_value=_post_resp(_experiment_data(), 201))
    get = MagicMock()
    monkeypatch.setattr(experiment_api.client, "post", post)
    monkeypatch.setattr(experiment_api.client, "get", get)

    before = datetime.now(timezone.utc)
    experiment, is_new = _call_create_or_resume()

    assert is_new is True
    assert post.call_args.args[0] == "/project/alice/demo/experiment"
    reported_created_at = post.call_args.args[1]["createdAt"]
    assert reported_created_at != "2024-08-01T00:00:00Z"
    # 新建实验直接复用上报的 createdAt，不做额外请求
    get.assert_not_called()
    assert experiment["createdAt"] == reported_created_at
    parsed = datetime.fromisoformat(reported_created_at.replace("Z", "+00:00"))
    assert before <= parsed <= datetime.now(timezone.utc)


def test_resume_skips_get_when_post_response_contains_created_at(monkeypatch):
    # 后端已按 SwanLab-Server#1090 在 POST 响应中携带 createdAt
    post = MagicMock(return_value=_post_resp(_experiment_data(createdAt="2024-07-01T00:00:00Z"), 200))
    get = MagicMock()
    monkeypatch.setattr(experiment_api.client, "post", post)
    monkeypatch.setattr(experiment_api.client, "get", get)

    experiment, is_new = _call_create_or_resume()

    assert is_new is False
    assert experiment["createdAt"] == "2024-07-01T00:00:00Z"
    get.assert_not_called()


def test_resume_falls_back_to_get_when_created_at_missing(monkeypatch):
    # 旧后端 POST 响应未携带 createdAt，回退到 GET 获取
    post = MagicMock(return_value=_post_resp(_experiment_data(), 200))
    get = MagicMock(return_value=SimpleNamespace(data=_experiment_data(createdAt="2024-07-01T00:00:00Z")))
    monkeypatch.setattr(experiment_api.client, "post", post)
    monkeypatch.setattr(experiment_api.client, "get", get)

    experiment, is_new = _call_create_or_resume()

    assert is_new is False
    assert experiment["createdAt"] == "2024-07-01T00:00:00Z"
    get.assert_called_once_with("/project/alice/demo/runs/experiment-id")


def test_resume_raises_when_get_also_missing_created_at(monkeypatch):
    post = MagicMock(return_value=_post_resp(_experiment_data(), 200))
    get = MagicMock(return_value=SimpleNamespace(data=_experiment_data()))
    monkeypatch.setattr(experiment_api.client, "post", post)
    monkeypatch.setattr(experiment_api.client, "get", get)

    with pytest.raises(ValueError, match="upgrade swanlab-server"):
        _call_create_or_resume()

    get.assert_called_once_with("/project/alice/demo/runs/experiment-id")


def _stale_timestamp() -> Timestamp:
    ts = Timestamp()
    ts.FromDatetime(datetime(2024, 8, 1, tzinfo=timezone.utc))
    return ts


def _stop_experiment(monkeypatch) -> MagicMock:
    put = MagicMock()
    monkeypatch.setattr(experiment_api.client, "put", put)
    return put


def test_stop_experiment_reports_current_time_as_finished_at(monkeypatch):
    # finishedAt 上报的是指标结束上报的时刻（当前时间），而非本地记录的结束时间
    put = _stop_experiment(monkeypatch)
    before = datetime.now(timezone.utc)

    experiment_api.stop_experiment(
        "alice", "demo", "experiment-id", state=RUN_STATE_FINISHED, finished_at=_stale_timestamp()
    )

    body = put.call_args.args[1]
    assert body["state"] == "FINISHED"
    assert body["from"] == "sdk"
    finished_at = datetime.fromisoformat(body["finishedAt"].replace("Z", "+00:00"))
    assert finished_at != datetime(2024, 8, 1, tzinfo=timezone.utc)
    assert before <= finished_at <= datetime.now(timezone.utc)


def test_stop_experiment_maps_crashed_and_aborted_states(monkeypatch):
    put = _stop_experiment(monkeypatch)

    experiment_api.stop_experiment(
        "alice", "demo", "experiment-id", state=RUN_STATE_ABORTED, finished_at=_stale_timestamp()
    )

    assert put.call_args.args[1]["state"] == "ABORTED"
