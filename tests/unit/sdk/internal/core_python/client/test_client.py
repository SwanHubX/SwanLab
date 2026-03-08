"""
@author: cunyue
@file: test_client.py
@time: 2026/3/7 20:45
@description: 测试 SwanLab 运行时客户端
"""

from datetime import datetime, timedelta, timezone
from unittest.mock import MagicMock, patch

import pytest
import responses

from swanlab.sdk.internal.core_python.client import (
    Client,
    _get_client,
    exists,
    new,
    reset,
)
from swanlab.sdk.internal.core_python.client import get as global_get
from swanlab.sdk.pkg.exceptions import ApiError


@pytest.fixture()
def mock_url():
    return "http://mock.example.com"


@pytest.fixture
def mock_login():
    patcher = patch("swanlab.sdk.internal.core_python.client.login_by_api_key")
    mock_func = patcher.start()
    mock_func.return_value = {"sid": "mock-token-123", "expiredAt": "2999-01-01T00:00:00.000Z"}
    yield mock_func
    patcher.stop()


@pytest.fixture
def client(mock_login, mock_url):
    _ = mock_login
    return Client(api_key="test-key", base_url=mock_url)


# -------------------------------------------------------------------
# Test Cases (测试用例)
# -------------------------------------------------------------------


def test_client_init_and_auth(mock_login, client, mock_url):
    """测试客户端初始化时，是否正确调用了登录接口并挂载 Cookie"""
    assert client._base_url == mock_url
    mock_login.assert_called_once_with(mock_url, "test-key", timeout=10)
    assert client._session.cookies.get("sid") == "mock-token-123"


def test_client_http_methods_and_url_join(client, mock_url):
    """测试 HTTP 请求的方法分发与 URL 拼接是否正确"""
    client._session.request = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = {"msg": "success"}
    client._session.request.return_value = mock_response

    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)

    # 1. 测试 GET 请求
    resp = client.get("/project", params={"id": 1})
    client._session.request.assert_called_with("GET", mock_url + "/project", params={"id": 1}, retries=None)
    assert resp.data == {"msg": "success"}  # 顺便验证数据类的 data 是否正确包装

    # 2. 测试 POST 请求
    client.post("run", data={"name": "test"})
    client._session.request.assert_called_with("POST", mock_url + "/run", json={"name": "test"}, retries=None)


def test_token_refresh_logic(client, mock_login):
    """测试当 token 即将过期时，发起请求是否会自动触发鉴权刷新"""
    assert mock_login.call_count == 1
    client._expired_at = datetime.now(timezone.utc)
    client._session.request = MagicMock()

    # 模拟正常响应防止 decode 报错
    mock_resp = MagicMock()
    mock_resp.json.return_value = {}
    client._session.request.return_value = mock_resp

    client.get("/test-refresh")
    assert mock_login.call_count == 2


def test_global_proxy_functions(mock_login, mock_url):
    """测试 new, exists, reset 全局状态的生命周期"""
    _ = mock_login

    # 清理初始状态以防干扰
    if exists():
        reset()

    # 1. 未初始化时调用，应该报错
    with pytest.raises(RuntimeError, match="not initialized"):
        _get_client()

    # 2. 正确初始化
    global_client = new("global-key", mock_url)
    assert exists() is True
    assert global_client._api_key == "global-key"

    # 3. 重复初始化应该报错
    with pytest.raises(RuntimeError, match="already exists"):
        new("another-key", mock_url)

    # 4. 测试全局快捷方法 (如 swanlab.client.get)
    global_client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    global_client._session.request = MagicMock()
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"ok": True}
    global_client._session.request.return_value = mock_resp

    resp = global_get("/global-test")
    global_client._session.request.assert_called_with("GET", mock_url + "/global-test", params=None, retries=None)
    assert resp.data == {"ok": True}

    # 5. 销毁并验证
    reset()
    assert exists() is False
    with pytest.raises(RuntimeError, match="not initialized"):
        reset()  # 再次销毁应报错


# -------------------------------------------------------------------
# Retry Mechanism Tests (重试机制测试)
# -------------------------------------------------------------------


@responses.activate()
def test_retry_default_on_server_error(client, mock_url):
    """默认重试：服务端持续返回 500，最终应抛出异常（而非静默失败）"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    responses.add(responses.GET, mock_url + "/health", status=500)

    # 优化：精确捕获 ApiError
    with pytest.raises(ApiError):
        client.get("/health")

    assert len(responses.calls) == 6


@responses.activate()
def test_retry_custom_zero_disables_retry(client, mock_url):
    """retries=0：禁用重试，第一次失败后立即抛出，调用次数恰好为 1"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    responses.add(responses.POST, mock_url + "/run", status=503)

    with pytest.raises(ApiError):
        client.post("/run", data={"name": "test"}, retries=0)

    assert len(responses.calls) == 1


@responses.activate()
def test_retry_custom_count(client, mock_url):
    """retries=2：前两次返回 500，第三次成功，最终应正常返回"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    target_url = mock_url + "/data"
    responses.add(responses.GET, target_url, status=500)
    responses.add(responses.GET, target_url, status=500)
    responses.add(responses.GET, target_url, json={"result": "ok"}, status=200)

    resp = client.get("/data", retries=2)

    assert resp.raw.status_code == 200
    assert resp.data == {"result": "ok"}
    assert len(responses.calls) == 3


def test_retry_invalid_negative_raises(client):
    """retries 为负数时，Adapter 应抛出 ValueError"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)

    with pytest.raises(ValueError, match="Invalid retry count"):
        client.get("/bad", retries=-1)


@responses.activate()
def test_retry_context_isolation(client, mock_url):
    """ContextVar 隔离：一次带 retries 的请求结束后，不影响下一次普通请求"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)

    responses.add(responses.GET, mock_url + "/a", status=500)
    responses.add(responses.GET, mock_url + "/b", json={"ok": True})

    with pytest.raises(ApiError):
        client.get("/a", retries=0)

    resp = client.get("/b")
    assert resp.raw.status_code == 200
    assert resp.data == {"ok": True}
