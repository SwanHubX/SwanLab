"""
@author: likang
@file: test_client_bootstrap.py
@time: 2026/11/3
@description: 登录握手错误分类与启动重试
"""

from unittest.mock import MagicMock

import pytest
import requests
import responses

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.sdk.internal.pkg.client import (
    is_login_auth_error,
    is_recoverable_login_error,
    is_transient_login_error,
)
from swanlab.sdk.internal.pkg.client.bootstrap import login_by_api_key, verify_api_key

BASE_URL = "http://api.test/api"
LOGIN_URL = f"{BASE_URL}/login/api_key"


def _wrap(cause: BaseException) -> AuthenticationError:
    try:
        raise AuthenticationError("wrapped") from cause
    except AuthenticationError as exc:
        return exc


@responses.activate
def test_login_propagates_http_error():
    """HTTP 错误不再被吞成 None。"""
    responses.add(responses.POST, LOGIN_URL, status=503)

    with pytest.raises(ApiError):
        login_by_api_key(BASE_URL, "test-key", timeout=1)

    # total=2，即首次请求 + 2 次重试
    assert len(responses.calls) == 3


@responses.activate
def test_login_propagates_network_error():
    """网络错误应保留原始异常类型。"""
    responses.add(responses.POST, LOGIN_URL, body=requests.ConnectionError("boom"))

    with pytest.raises(requests.ConnectionError):
        login_by_api_key(BASE_URL, "test-key", timeout=1)


@responses.activate
def test_login_retries_transient_5xx_then_succeeds():
    """瞬时 5xx 允许重试后成功。"""
    responses.add(responses.POST, LOGIN_URL, status=500)
    responses.add(responses.POST, LOGIN_URL, json={"sid": "mock-sid"}, status=200)

    assert login_by_api_key(BASE_URL, "test-key", timeout=1) == {"sid": "mock-sid"}
    assert len(responses.calls) == 2


@responses.activate
def test_verify_api_key_returns_response():
    """校验成功时返回原始登录响应。"""
    responses.add(responses.POST, LOGIN_URL, json={"sid": "mock-sid"}, status=200)

    assert verify_api_key(BASE_URL, "test-key", timeout=1) == {"sid": "mock-sid"}


@responses.activate
def test_verify_api_key_returns_none_on_error():
    """校验失败时返回 None，供 CLI 展示统一提示。"""
    responses.add(responses.POST, LOGIN_URL, status=503)

    assert verify_api_key(BASE_URL, "test-key", timeout=1) is None


def test_transient_login_error_unwraps_client_cause():
    """Client 包装的一层 __cause__ 应被识别为可降级的瞬时错误。"""
    assert is_transient_login_error(_wrap(requests.ConnectionError("boom")))


def test_unauthorized_login_error_is_not_transient():
    """401 应识别为认证错误，可提示用户重新输入 API Key。"""
    cause = ApiError(MagicMock(status_code=401), method="POST", code=401, message="unauthorized")
    exc = _wrap(cause)

    assert not is_transient_login_error(exc)
    assert is_login_auth_error(exc)


@responses.activate
def test_non_json_login_response_is_transient_but_not_recoverable():
    """2xx 但不是 JSON（如强制门户）不能进入后台恢复流程。"""
    responses.add(
        responses.POST,
        LOGIN_URL,
        body="<html>captive portal</html>",
        status=200,
        content_type="text/html",
    )

    with pytest.raises(requests.JSONDecodeError) as exc_info:
        login_by_api_key(BASE_URL, "test-key", timeout=1)
    exc = _wrap(exc_info.value)

    assert is_transient_login_error(exc)
    assert not is_recoverable_login_error(exc)


def test_invalid_login_payload_is_transient_but_not_recoverable():
    """2xx 但缺少 sid/expiredAt/用户信息属于协议错误，不能进入后台恢复。"""
    cause = ApiError(None, method="POST", code="invalid_login_response", message="missing sid", retryable=True)
    exc = _wrap(cause)

    assert is_transient_login_error(exc)
    assert not is_recoverable_login_error(exc)


def test_recoverable_login_error_keeps_network_and_5xx():
    """网络错误、408/429/5xx 仍允许后台恢复；401/403 与配置错误不允许。"""
    assert is_recoverable_login_error(_wrap(requests.ConnectionError("boom")))
    for status in (408, 429, 500, 503):
        assert is_recoverable_login_error(_wrap(ApiError(MagicMock(status_code=status), method="POST", code=status)))
    for status in (401, 403):
        assert not is_recoverable_login_error(
            _wrap(ApiError(MagicMock(status_code=status), method="POST", code=status))
        )
    assert not is_recoverable_login_error(_wrap(requests.exceptions.InvalidURL("bad url")))


def test_login_error_helpers_tolerate_cyclic_cause():
    """即使异常链被手工构造成环，也只解包一层，不会死循环。"""
    exc = AuthenticationError("cyclic")
    exc.__cause__ = exc

    assert not is_transient_login_error(exc)
    assert not is_login_auth_error(exc)
