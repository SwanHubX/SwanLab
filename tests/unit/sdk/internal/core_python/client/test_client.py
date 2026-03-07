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
    get_client,
    init_client,
    reset_client,
)
from swanlab.sdk.internal.core_python.client import (
    get as global_get,
)


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
    # 1. 验证是否正确去除了 base_url 末尾的斜杠
    assert client._base_url == mock_url

    # 2. 验证是否在初始化时就调用了登录逻辑
    mock_login.assert_called_once_with(mock_url, "test-key")

    # 3. 验证 session 中是否成功注入了 sid
    assert client._session.cookies.get("sid") == "mock-token-123"


def test_client_http_methods_and_url_join(client, mock_url):
    """测试 HTTP 请求的方法分发与 URL 拼接是否正确"""
    # Mock 掉底层的 session.request，防止发出真实网络请求
    client._session.request = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = {"msg": "success"}
    client._session.request.return_value = mock_response

    # 给个安全的未过期时间，防止触发 _refresh_auth
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)

    # 1. 测试 GET 请求，URL 带前导斜杠
    client.get("/project", params={"id": 1})
    client._session.request.assert_called_with("GET", mock_url + "/project", params={"id": 1}, retries=None)

    # 2. 测试 POST 请求，URL 不带前导斜杠（验证 lstrip 的效果）
    client.post("run", data={"name": "test"})
    client._session.request.assert_called_with("POST", mock_url + "/run", json={"name": "test"}, retries=None)


def test_token_refresh_logic(client, mock_login):
    """测试当 token 即将过期时，发起请求是否会自动触发鉴权刷新"""
    # 初始化时已经调过一次了
    assert mock_login.call_count == 1

    # 强制将过期时间设置为【当前时间】（意味着已经处于 REFRESH_TIME 缓冲期内）
    client._expired_at = datetime.now(timezone.utc)

    # Mock 发送请求
    client._session.request = MagicMock()
    client._session.request.return_value = MagicMock()

    # 触发随便一个请求
    client.get("/test-refresh")

    # 验证是否再次调用了登录接口进行刷新
    assert mock_login.call_count == 2


def test_global_proxy_functions(mock_login, mock_url):
    """测试 init_client, get_client, reset_client 全局状态的生命周期"""
    _ = mock_login
    # 1. 未初始化时调用，应该报错
    reset_client()
    with pytest.raises(RuntimeError, match="not initialized"):
        get_client()

    # 2. 正确初始化
    init_client("global-key", mock_url)
    global_client = get_client()
    assert global_client._api_key == "global-key"

    # 3. 测试全局快捷方法 (如 swanlab.client.get)
    global_client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    global_client._session.request = MagicMock()
    global_client._session.request.return_value = MagicMock()

    global_get("/global-test")
    global_client._session.request.assert_called_with("GET", mock_url + "/global-test", params=None, retries=None)


# -------------------------------------------------------------------
# Retry Mechanism Tests (重试机制测试)
# -------------------------------------------------------------------


@responses.activate()
def test_retry_default_on_server_error(client, mock_url):
    """默认重试：服务端持续返回 500，最终应抛出异常（而非静默失败）"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    # 注册一个始终返回 500 的端点
    responses.add(responses.GET, mock_url + "/health", status=500)

    with pytest.raises(Exception):
        client.get("/health")

    assert len(responses.calls) == 6


@responses.activate()
def test_retry_custom_zero_disables_retry(client, mock_url):
    """retries=0：禁用重试，第一次失败后立即抛出，调用次数恰好为 1"""
    client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    responses.add(responses.POST, mock_url + "/run", status=503)

    with pytest.raises(Exception):
        client.post("/run", data={"name": "test"}, retries=0)

    # 优化：验证 calls 的长度
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

    # 验证总共确实发生了 3 次网络调用
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

    # 第一次请求带 retries=0，失败抛出
    with pytest.raises(Exception):
        client.get("/a", retries=0)

    # 第二次请求不带 retries，ContextVar 应已被清理，正常走默认重试策略并成功
    resp = client.get("/b")
    assert resp.raw.status_code == 200
    assert resp.data == {"ok": True}
