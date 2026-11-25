"""
@author: Zirui Cai
@file: test_checkpoint_e2e.py
@time: 2026/11/25
@description: swanlab.checkpoint() 端到端测试
"""

import json
from typing import cast

import pytest
import responses as responses_lib

from swanlab.sdk.cmd.checkpoint import checkpoint
from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.login import login_cli
from swanlab.sdk.cmd.merge_settings import merge_settings
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import get_run, has_run
from swanlab.sdk.internal.settings import create_settings

# ============================================================
# 常量
# ============================================================

API_HOST = "https://api.fake.swanlab.cn"
WEB_HOST = "https://test.swanlab.cn"
USERNAME = "test-user"
PROJECT = "test-project"
RUN_ID = "test-run-id"
EXPERIMENT_CUID = "test-experiment-cuid"
API_KEY = "test-api-key"
CHECKPOINT_CUID = "test-checkpoint-cuid"
URI = "s3://my-bucket/run-a/checkpoint-1000"


# ============================================================
# 响应体工厂
# ============================================================


def make_login_resp(**overrides) -> dict:
    return {
        "sid": "mock-sid",
        "expiredAt": "2099-12-31T23:59:59.000Z",
        "userInfo": {"username": USERNAME, "name": "Test User"},
        **overrides,
    }


def make_project_detail_resp(**overrides) -> dict:
    return {
        "cuid": "test-project-cuid",
        "name": PROJECT,
        "version": 1,
        "group": {"username": USERNAME},
        "username": USERNAME,
        "path": f"/{USERNAME}/{PROJECT}",
        "visibility": "PRIVATE",
        "_count": {"experiments": 0, "contributors": 1, "collaborators": 0, "clones": 0},
        **overrides,
    }


def make_experiment_resp(**overrides) -> dict:
    return {"cuid": EXPERIMENT_CUID, "slug": RUN_ID, "name": "test-experiment", **overrides}


def make_ok_resp(**overrides) -> dict:
    return {"message": "ok", **overrides}


# ============================================================
# HTTP Mock Fixtures
# ============================================================


@pytest.fixture
def rsps():
    with responses_lib.RequestsMock() as r:
        yield r
        if has_run():
            get_run().finish()


@pytest.fixture
def mock_online_settings():
    merge_settings({"api_host": API_HOST, "web_host": WEB_HOST, "probe": {"monitor": False}})


@pytest.fixture
def mock_login_api(rsps):
    rsps.add(responses_lib.POST, f"{API_HOST}/api/login/api_key", json=make_login_resp(), status=200)
    return rsps


@pytest.fixture
def mock_project_get_api(rsps):
    rsps.add(
        responses_lib.GET,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}",
        json=make_project_detail_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_experiment_create_api(rsps):
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/experiment",
        json=make_experiment_resp(),
        status=201,
    )
    return rsps


@pytest.fixture
def mock_experiment_stop_api(rsps):
    rsps.add(
        responses_lib.PUT,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/runs/{EXPERIMENT_CUID}/state",
        json=make_ok_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_profile_api(rsps):
    rsps.add(
        responses_lib.PUT,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/runs/{EXPERIMENT_CUID}/profile",
        json=make_ok_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_columns_api(rsps):
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/projects/{USERNAME}/{PROJECT}/series",
        json=make_ok_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_metrics_api(rsps):
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/house/metrics",
        json=make_ok_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_heartbeat_api(rsps):
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/house/experiments/{EXPERIMENT_CUID}/heartbeat",
        json=make_ok_resp(),
        status=200,
    )
    return rsps


@pytest.fixture
def mock_checkpoint_api(rsps):
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/checkpoint",
        json={"cuid": CHECKPOINT_CUID},
        status=201,
    )
    return rsps


@pytest.fixture
def logged_in_client(mock_login_api, mock_online_settings):
    login_cli(api_key=API_KEY, host=create_settings().api_host)
    yield


@pytest.fixture
def mock_online_init_apis(
    logged_in_client,
    mock_project_get_api,
    mock_experiment_create_api,
    mock_experiment_stop_api,
    mock_profile_api,
    mock_metrics_api,
    mock_heartbeat_api,
):
    pass


@pytest.fixture
def mock_online_checkpoint_apis(mock_online_init_apis, mock_checkpoint_api):
    pass


@pytest.fixture
def mock_online_checkpoint_apis_with_log(
    mock_online_init_apis, mock_columns_api, mock_checkpoint_api
):
    """run.log() を含むテスト用：series エンドポイントも登録する"""
    pass


# ============================================================
# 非 online 模式：直接用 mode= 参数，不需要任何 HTTP mock
# ============================================================


class TestCheckpointNonOnlineMode:
    def test_disabled_returns_empty_and_warns(self, monkeypatch):
        """disabled 模式：warning + 返回空字符串"""
        warnings = []
        monkeypatch.setattr(console, "warning", lambda *a, **kw: warnings.append(str(a)))
        with init(mode="disabled"):
            result = checkpoint(uri=URI)
        assert result == ""
        assert any("online mode" in w for w in warnings)

    def test_local_returns_empty_and_warns(self, monkeypatch):
        """local 模式：warning + 返回空字符串"""
        warnings = []
        monkeypatch.setattr(console, "warning", lambda *a, **kw: warnings.append(str(a)))
        with init(mode="local"):
            result = checkpoint(uri=URI)
        assert result == ""
        assert any("online mode" in w for w in warnings)

    def test_offline_returns_empty(self):
        """offline 模式：返回空字符串"""
        with init(mode="offline"):
            assert checkpoint(uri=URI) == ""

    def test_no_run_no_client_returns_empty(self, monkeypatch):
        """无活跃 run 且无 client：返回空字符串"""
        warnings = []
        monkeypatch.setattr(console, "warning", lambda *a, **kw: warnings.append(str(a)))
        result = checkpoint(uri=URI)
        assert result == ""
        assert any("run" in w for w in warnings)

    def test_does_not_raise_in_any_non_online_mode(self):
        """任何非 online 场景下不抛出异常"""
        with init(mode="disabled"):
            try:
                checkpoint(uri=URI, step=None, labels=["best"], metadata={"lr": 0.01})
            except Exception as e:
                pytest.fail(f"checkpoint() raised: {e}")


# ============================================================
# Online 模式
# ============================================================


class TestCheckpointOnlineMode:
    def test_returns_cuid_on_success(self, mock_online_checkpoint_apis):
        """正常登记返回平台分配的 cuid"""
        init(mode="online", project=PROJECT)
        assert checkpoint(uri=URI) == CHECKPOINT_CUID

    def test_payload_contains_experiment_id_and_uri(self, mock_online_checkpoint_apis, rsps):
        """请求 body 包含 experimentId 和 uri"""
        init(mode="online", project=PROJECT)
        checkpoint(uri=URI, name="best", step=1000)

        calls = [c for c in rsps.calls if c.request.url.endswith("/api/checkpoint")]
        assert len(calls) == 1
        body = json.loads(cast(bytes, calls[0].request.body))
        assert body["experimentId"] == EXPERIMENT_CUID
        assert body["uri"] == URI
        assert body["name"] == "best"
        assert body["step"] == 1000

    def test_metadata_and_labels_forwarded(self, mock_online_checkpoint_apis, rsps):
        """metadata 和 labels 正确传递"""
        init(mode="online", project=PROJECT)
        checkpoint(uri=URI, metadata={"optimizer": "AdamW"}, labels=["best", "ready"])

        body = json.loads(cast(bytes, [
            c for c in rsps.calls if c.request.url.endswith("/api/checkpoint")
        ][0].request.body))
        assert body["metadata"] == {"optimizer": "AdamW"}
        assert body["labels"] == [{"name": "best"}, {"name": "ready"}]

    def test_server_error_returns_empty_string(self, mock_online_init_apis, rsps):
        """后端 4xx 返回空字符串，不抛异常"""
        rsps.add(
            responses_lib.POST,
            f"{API_HOST}/api/checkpoint",
            json={"message": "bad request"},
            status=400,
        )
        init(mode="online", project=PROJECT)
        assert checkpoint(uri=URI) == ""


# ============================================================
# Step 解析
# ============================================================


class TestCheckpointStepResolution:
    def _get_checkpoint_body(self, rsps) -> dict:
        calls = [c for c in rsps.calls if c.request.url.endswith("/api/checkpoint")]
        return json.loads(cast(bytes, calls[0].request.body))

    def test_auto_step_from_active_run(self, mock_online_checkpoint_apis_with_log, rsps):
        """未传 step：从 active run 的 global_step 取"""
        run = init(mode="online", project=PROJECT)
        run.log({"loss": 0.5})   # step 0
        run.log({"loss": 0.4})   # step 1
        checkpoint(uri=URI)
        assert self._get_checkpoint_body(rsps)["step"] == 1

    def test_explicit_step_overrides_auto(self, mock_online_checkpoint_apis_with_log, rsps):
        """显式传 step 时使用传入值"""
        run = init(mode="online", project=PROJECT)
        run.log({"loss": 0.5})
        checkpoint(uri=URI, step=500)
        assert self._get_checkpoint_body(rsps)["step"] == 500

    def test_explicit_none_step_sends_null(self, mock_online_checkpoint_apis_with_log, rsps):
        """step=None 时 body 中为 null"""
        run = init(mode="online", project=PROJECT)
        run.log({"loss": 0.5})
        checkpoint(uri=URI, step=None)
        assert self._get_checkpoint_body(rsps).get("step") is None

    def test_no_log_step_is_none(self, mock_online_checkpoint_apis, rsps):
        """run 从未 log（global_step=-1）时 step 为 null"""
        init(mode="online", project=PROJECT)
        checkpoint(uri=URI)
        assert self._get_checkpoint_body(rsps).get("step") is None
