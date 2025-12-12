"""
@author: cunyue
@file: test_session.py
@time: 2025/9/9 15:12
@description: $END$
"""

from unittest.mock import Mock, patch

import pytest
import requests
import responses
from requests.adapters import HTTPAdapter
from responses import registries
from urllib3.util.retry import Retry

from swanlab.core_python import create_session
from swanlab.core_python.client.session import TimeoutHTTPAdapter, DEFAULT_TIMEOUT
from swanlab.package import get_package_version


@pytest.mark.parametrize("url", ["https://api.example.com/retry", "http://api.example.com/retry"])
@responses.activate(registry=registries.OrderedRegistry)
def test_retry(url):
    """
    测试重试机制
    """

    [responses.add(responses.GET, url, body="Error", status=500) for _ in range(5)]
    responses.add(responses.GET, url, body="Success", status=200)
    s = create_session()
    resp = s.get(url)
    assert resp.text == "Success"
    assert len(responses.calls) == 6


@responses.activate(registry=registries.OrderedRegistry)
def test_session_headers():
    """
    测试会话是否包含正确的自定义请求头
    """
    # 1. 准备测试数据
    test_url = "https://api.example.com/test"
    expected_sdk_version = get_package_version()

    # 2. 模拟响应 - 捕获请求头
    captured_headers = {}

    def request_callback(request):
        # 捕获所有请求头
        nonlocal captured_headers
        captured_headers = dict(request.headers)
        return (200, {}, "OK")

    responses.add_callback(responses.GET, test_url, callback=request_callback)

    # 3. 创建会话并发送请求
    session = create_session()
    response = session.get(test_url)

    # 4. 验证
    assert response.status_code == 200

    # 验证自定义头存在且值正确
    assert "swanlab-sdk" in captured_headers
    assert captured_headers["swanlab-sdk"] == expected_sdk_version

    # 验证User-Agent等默认头也存在（可选）
    assert "User-Agent" in captured_headers

    # 打印所有捕获的请求头（调试用）
    print("\n捕获的请求头:", captured_headers)


@responses.activate(registry=registries.OrderedRegistry)
def test_header_merging():
    """
    测试请求级别headers与会话级别headers的合并
    """
    test_url = "https://api.example.com/merge"
    custom_header = {"X-Custom-Request-Header": "test-value"}

    captured_headers = {}

    def request_callback(request):
        nonlocal captured_headers
        captured_headers = dict(request.headers)
        return (200, {}, "OK")

    responses.add_callback(responses.GET, test_url, callback=request_callback)

    # 创建会话（自带swanlab-sdk头）
    session = create_session()

    # 发送带额外请求头的请求
    response = session.get(test_url, headers=custom_header)

    # 验证
    assert response.status_code == 200

    # 验证会话头依然存在
    assert "swanlab-sdk" in captured_headers

    # 验证请求级别头已添加
    assert "X-Custom-Request-Header" in captured_headers
    assert captured_headers["X-Custom-Request-Header"] == "test-value"

    # 验证合并而非覆盖（两个头都存在）
    assert len(captured_headers) >= 2


# ========== TimeoutHTTPAdapter Tests ==========


def test_timeout_adapter_initialization_with_timeout():
    """
    测试TimeoutHTTPAdapter初始化时正确设置timeout参数
    """
    timeout_value = 30
    adapter = TimeoutHTTPAdapter(timeout=timeout_value)

    assert adapter.timeout == timeout_value


def test_timeout_adapter_initialization_without_timeout():
    """
    测试TimeoutHTTPAdapter初始化时未提供timeout参数的情况
    """
    adapter = TimeoutHTTPAdapter()

    assert adapter.timeout is None


def test_timeout_adapter_initialization_with_other_params():
    """
    测试TimeoutHTTPAdapter初始化时传递其他HTTPAdapter参数
    """
    retry = Retry(total=3)
    adapter = TimeoutHTTPAdapter(max_retries=retry, timeout=45)

    assert adapter.timeout == 45
    assert adapter.max_retries == retry


def test_timeout_adapter_uses_default_timeout():
    """
    测试TimeoutHTTPAdapter在未显式指定timeout时使用默认timeout
    """
    test_url = "https://api.example.com/timeout-test"

    # 创建adapter并挂载到session
    adapter = TimeoutHTTPAdapter(timeout=25)

    session = requests.Session()
    session.mount("https://", adapter)

    # 创建一个PreparedRequest来测试send方法
    req = requests.Request('GET', test_url)
    prepared = session.prepare_request(req)

    # 直接调用adapter的send方法，不传timeout
    # 我们期望adapter会自动添加timeout=25
    with patch.object(HTTPAdapter, 'send') as mock_parent_send:
        mock_parent_send.return_value = Mock(status_code=200, text="OK")

        adapter.send(prepared)

        # 验证父类的send方法被调用时包含了timeout参数
        call_kwargs = mock_parent_send.call_args[1]
        assert 'timeout' in call_kwargs
        assert call_kwargs['timeout'] == 25


def test_timeout_adapter_respects_explicit_timeout():
    """
    测试TimeoutHTTPAdapter在显式指定timeout时覆盖默认timeout
    """
    test_url = "https://api.example.com/explicit-timeout"

    # 创建adapter，设置默认timeout为30
    adapter = TimeoutHTTPAdapter(timeout=30)

    session = requests.Session()
    session.mount("https://", adapter)

    # 创建一个PreparedRequest来测试send方法
    req = requests.Request('GET', test_url)
    prepared = session.prepare_request(req)

    # 直接调用adapter的send方法，显式传timeout=10
    with patch.object(HTTPAdapter, 'send') as mock_parent_send:
        mock_parent_send.return_value = Mock(status_code=200, text="OK")

        adapter.send(prepared, timeout=10)

        # 验证父类的send方法被调用时使用了显式指定的timeout=10
        call_kwargs = mock_parent_send.call_args[1]
        assert 'timeout' in call_kwargs
        assert call_kwargs['timeout'] == 10


def test_timeout_adapter_with_none_timeout():
    """
    测试TimeoutHTTPAdapter当timeout为None时不注入超时
    """
    test_url = "https://api.example.com/none-timeout"

    # 创建adapter，timeout为None
    adapter = TimeoutHTTPAdapter(timeout=None)

    session = requests.Session()
    session.mount("https://", adapter)

    # 创建一个PreparedRequest来测试send方法
    req = requests.Request('GET', test_url)
    prepared = session.prepare_request(req)

    # 直接调用adapter的send方法，不传timeout
    with patch.object(HTTPAdapter, 'send') as mock_parent_send:
        mock_parent_send.return_value = Mock(status_code=200, text="OK")

        adapter.send(prepared)

        # 验证父类的send方法被调用时不应包含timeout参数
        call_kwargs = mock_parent_send.call_args[1]
        assert 'timeout' not in call_kwargs


@responses.activate
def test_create_session_uses_default_timeout():
    """
    测试create_session创建的会话使用DEFAULT_TIMEOUT
    """
    test_url = "https://api.example.com/session-timeout"

    responses.add(responses.GET, test_url, body="OK", status=200)

    session = create_session()

    # 获取adapter
    adapter = session.get_adapter(test_url)
    assert isinstance(adapter, TimeoutHTTPAdapter)
    assert adapter.timeout == DEFAULT_TIMEOUT

    # 创建一个PreparedRequest来测试send方法
    req = requests.Request('GET', test_url)
    prepared = session.prepare_request(req)

    # 验证实际请求使用该timeout
    with patch.object(HTTPAdapter, 'send') as mock_parent_send:
        mock_parent_send.return_value = Mock(status_code=200, text="OK")

        adapter.send(prepared)

        # 验证使用了DEFAULT_TIMEOUT
        call_kwargs = mock_parent_send.call_args[1]
        assert 'timeout' in call_kwargs
        assert call_kwargs['timeout'] == DEFAULT_TIMEOUT


def test_timeout_adapter_inherits_from_httpAdapter():
    """
    测试TimeoutHTTPAdapter正确继承自HTTPAdapter
    """
    adapter = TimeoutHTTPAdapter(timeout=20)

    assert isinstance(adapter, HTTPAdapter)
    assert isinstance(adapter, TimeoutHTTPAdapter)
