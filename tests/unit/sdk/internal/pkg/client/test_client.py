"""
@author: cunyue
@file: test_client.py
@time: 2026/4/14 00:44
@description: 测试 SwanLab 运行时客户端的核心功能，包括初始化、鉴权、请求发送和重试。
"""

from unittest.mock import MagicMock

import pytest
import requests
import responses

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.sdk.internal.pkg.client import Client, session, verify_api_key

PROFILE = {
    "uid": 1,
    "avatar": "",
    "username": "test-user",
    "name": "Test User",
    "createdAt": "2026-01-01T00:00:00.000Z",
    "verified": True,
}


@pytest.fixture()
def mock_base_url():
    """用户传入的基础 URL"""
    return "http://mock.example.com"


@pytest.fixture()
def mock_api_url(mock_base_url):
    """Client 内部强制拼接后的真实 API URL"""
    return f"{mock_base_url}/api"


@pytest.fixture()
def authenticated_client(mock_base_url, mock_api_url):
    with responses.RequestsMock() as rsps:
        rsps.add(responses.GET, f"{mock_api_url}/auth/verify", json=PROFILE)
        client = Client(api_key="test-key", base_url=mock_base_url)
        try:
            yield client, rsps
        finally:
            client.close()


@pytest.fixture()
def session_close_spy(monkeypatch):
    http = session.create()
    close = MagicMock(wraps=http.close)
    monkeypatch.setattr(http, "close", close)
    monkeypatch.setattr(session, "create", lambda: http)
    try:
        yield close
    finally:
        http.close()


def test_client_init_and_auth(authenticated_client, mock_api_url):
    """初始化时挂载 ApiKey 常驻头，并通过 /verify 缓存用户档案"""
    client, _ = authenticated_client

    assert client._base_url == mock_api_url
    assert client._session.headers["Authorization"] == "ApiKey test-key"
    assert client.profile["username"] == "test-user"
    assert client.username == "test-user"


def test_client_http_methods_and_url_join(authenticated_client, mock_api_url):
    """测试 HTTP 请求的方法分发与 URL 拼接是否正确"""
    client, _ = authenticated_client
    client._session.request = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = {"msg": "success"}
    client._session.request.return_value = mock_response

    # 1. 测试 GET 请求
    resp = client.get("/project", params={"id": 1})
    client._session.request.assert_called_with(
        "GET", mock_api_url + "/project", params={"id": 1}, retries=None, log_error=True
    )
    assert resp.data == {"msg": "success"}  # 顺便验证数据类的 data 是否正确包装

    # 2. 测试 POST 请求
    client.post("run", data={"name": "test"})
    client._session.request.assert_called_with(
        "POST", mock_api_url + "/run", json={"name": "test"}, retries=None, log_error=True
    )

    # 3. 测试单次请求错误日志开关透传
    client.post("run", data={"name": "test"}, log_error=False)
    client._session.request.assert_called_with(
        "POST", mock_api_url + "/run", json={"name": "test"}, retries=None, log_error=False
    )


# -------------------------------------------------------------------
# Retry Mechanism Tests (重试机制测试)
# -------------------------------------------------------------------


def test_retry_default_on_server_error(authenticated_client, mock_api_url):
    """默认重试：服务端持续返回 500，最终应抛出异常（而非静默失败）"""
    client, rsps = authenticated_client

    rsps.add(responses.GET, mock_api_url + "/health", status=500)

    with pytest.raises(ApiError):
        client.get("/health")

    # verify(1) + health(1 + 5 次默认重试)
    assert len(rsps.calls) == 7


def test_retry_custom_zero_disables_retry(authenticated_client, mock_api_url):
    """retries=0：禁用重试，第一次失败后立即抛出"""
    client, rsps = authenticated_client

    rsps.add(responses.POST, mock_api_url + "/run", status=503)

    with pytest.raises(ApiError):
        client.post("/run", data={"name": "test"}, retries=0)

    # verify(1) + run(1)
    assert len(rsps.calls) == 2


def test_retry_custom_count(authenticated_client, mock_api_url):
    """retries=2：前两次返回 500，第三次成功，最终应正常返回"""
    client, rsps = authenticated_client
    target_url = mock_api_url + "/data"

    rsps.add(responses.GET, target_url, status=500)
    rsps.add(responses.GET, target_url, status=500)
    rsps.add(responses.GET, target_url, json={"result": "ok"}, status=200)

    resp = client.get("/data", retries=2)

    assert resp.raw.status_code == 200
    assert resp.data == {"result": "ok"}
    # verify(1) + data(3)
    assert len(rsps.calls) == 4


def test_retry_invalid_negative_raises(authenticated_client, mock_api_url):
    """retries 为负数时，Adapter 应抛出 ValueError"""
    client, _ = authenticated_client

    with pytest.raises(ValueError, match="Invalid retry count"):
        client.get("/bad", retries=-1)


def test_retry_context_isolation(authenticated_client, mock_api_url):
    """ContextVar 隔离：一次带 retries 的请求结束后，不影响下一次普通请求"""
    client, rsps = authenticated_client

    rsps.add(responses.GET, mock_api_url + "/a", status=500)
    rsps.add(responses.GET, mock_api_url + "/b", json={"ok": True})

    with pytest.raises(ApiError):
        client.get("/a", retries=0)

    resp = client.get("/b")
    assert resp.raw.status_code == 200
    assert resp.data == {"ok": True}
    # verify(1) + a(1) + b(1)
    assert len(rsps.calls) == 3


@pytest.mark.parametrize("status", [401, 403, 404, 429, 503])
@responses.activate
def test_verification_error_preserves_category_and_closes_session(
    status, session_close_spy, mock_base_url, mock_api_url
):
    responses.add(
        responses.GET,
        f"{mock_api_url}/auth/verify",
        status=status,
        json={"message": "user is not verified"},
        headers={"Retry-After": "5"},
    )
    error_type = AuthenticationError if status in (401, 403) else RuntimeError if status == 404 else ApiError
    with pytest.raises(error_type) as caught:
        verify_api_key(mock_base_url, "test-key")
    session_close_spy.assert_called_once()
    assert len(responses.calls) == 1
    if status in (401, 403):
        assert "user is not verified" in str(caught.value)
    elif status == 404:
        assert "Please upgrade" in str(caught.value)
    else:
        assert isinstance(caught.value, ApiError)
        assert caught.value.response.status_code == status
        assert caught.value.response.headers["Retry-After"] == "5"


@pytest.mark.parametrize("body", ["", "<html>wrong route</html>"])
@responses.activate
def test_empty_unauthorized_response_has_fallback(body, mock_base_url, mock_api_url):
    responses.add(responses.GET, f"{mock_api_url}/auth/verify", status=401, body=body)
    with pytest.raises(AuthenticationError, match="invalid API key"):
        Client("test-key", mock_base_url)


@pytest.mark.parametrize(
    "data",
    [
        None,
        [],
        "<html>wrong route</html>",
        {},
        {**PROFILE, "uid": True},
        {**PROFILE, "uid": 0},
        {**PROFILE, "username": " "},
        {**PROFILE, "createdAt": None},
        {**PROFILE, "verified": "true"},
        {**PROFILE, "name": None},
    ],
)
@responses.activate
def test_malformed_profile_closes_session(data, session_close_spy, mock_base_url, mock_api_url):
    responses.add(responses.GET, f"{mock_api_url}/auth/verify", status=200, json=data)
    with pytest.raises(ValueError, match="Invalid authentication response"):
        Client("test-key", mock_base_url)
    session_close_spy.assert_called_once()


@responses.activate
def test_verify_accepts_omitted_optional_fields_and_closes(session_close_spy, mock_base_url, mock_api_url):
    responses.add(
        responses.GET,
        f"{mock_api_url}/auth/verify",
        json={
            "uid": 1,
            "username": "test-user",
            "createdAt": PROFILE["createdAt"],
        },
    )
    profile = verify_api_key(mock_base_url, "test-key")
    assert profile["name"] == ""
    assert profile["verified"] is False
    session_close_spy.assert_called_once()


@pytest.mark.parametrize("failure", [requests.Timeout("timeout"), KeyboardInterrupt()])
def test_verification_interruption_closes_session(failure, monkeypatch, mock_base_url):
    http = MagicMock()
    http.request.side_effect = failure
    monkeypatch.setattr("swanlab.sdk.internal.pkg.client.session.create", lambda: http)
    with pytest.raises(type(failure)):
        verify_api_key(mock_base_url, "test-key")
    http.close.assert_called_once()
