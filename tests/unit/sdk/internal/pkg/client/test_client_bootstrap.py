"""
@author: likang
@file: test_client_bootstrap.py
@time: 2026/11/3
@description: 登录握手错误分类与快速失败
"""

import pytest
import requests
import responses

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.pkg.client.bootstrap import login_by_api_key, verify_api_key

BASE_URL = "http://api.test/api"
LOGIN_URL = f"{BASE_URL}/login/api_key"


@responses.activate
def test_login_propagates_http_error():
    """HTTP 错误不再被吞成 None。"""
    responses.add(responses.POST, LOGIN_URL, status=503)

    with pytest.raises(ApiError):
        login_by_api_key(BASE_URL, "test-key", timeout=1)

    # total=1，即首次请求 + 1 次重试
    assert len(responses.calls) == 2


@responses.activate
def test_login_propagates_network_error():
    """网络错误应保留原始异常类型。"""
    responses.add(responses.POST, LOGIN_URL, body=requests.ConnectionError("boom"))

    with pytest.raises(requests.ConnectionError):
        login_by_api_key(BASE_URL, "test-key", timeout=1)


@responses.activate
def test_login_retries_once_then_succeeds():
    """瞬时 5xx 允许重试一次。"""
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
