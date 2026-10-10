"""
@author: cunyue
@file: test_core_client.py
@time: 2026/3/7 20:45
@description: 测试 SwanLab 运行时全局客户端单例生命周期
"""

import multiprocessing
from unittest.mock import MagicMock

import pytest
import responses

from swanlab.sdk.internal.core_python.client import (
    _get_client,
    exists,
    new,
    reset,
)
from swanlab.sdk.internal.core_python.client import get as global_get
from swanlab.sdk.internal.core_python.client import post as global_post
from swanlab.sdk.internal.core_python.client import username as global_username

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


# -------------------------------------------------------------------
# Test Cases (测试用例)
# -------------------------------------------------------------------


@responses.activate
def test_global_proxy_functions(mock_base_url, mock_api_url):
    """测试 new, exists, reset 全局状态的生命周期"""
    responses.add(responses.GET, mock_api_url + "/auth/verify", json=PROFILE, status=200)

    # 清理初始状态以防干扰
    if exists():
        reset()

    # 1. 未初始化时调用，应该报错
    with pytest.raises(RuntimeError, match="not initialized"):
        _get_client()

    # 2. 正确初始化
    global_client = new("global-key", mock_base_url)
    assert exists() is True
    assert responses.calls[0].request.headers["Authorization"] == "ApiKey global-key"
    assert global_username() == "test-user"

    # 3. 重复初始化应该报错
    with pytest.raises(RuntimeError, match="already exists"):
        new("another-key", mock_base_url)

    # 4. 测试全局快捷方法 (如 swanlab.client.get)
    global_client._session.request = MagicMock()
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"ok": True}
    global_client._session.request.return_value = mock_resp

    resp = global_get("/global-test")

    # 注意这里断言的应该是拼接后的 mock_api_url
    global_client._session.request.assert_called_with(
        "GET", mock_api_url + "/global-test", params=None, retries=None, log_error=True
    )
    assert resp.data == {"ok": True}

    global_post("/global-submit", data={"name": "test"}, log_error=False)
    global_client._session.request.assert_called_with(
        "POST", mock_api_url + "/global-submit", json={"name": "test"}, retries=None, log_error=False
    )

    # 5. 销毁并验证
    close = MagicMock(wraps=global_client.close)
    global_client.close = close
    reset()
    close.assert_called_once()
    assert exists() is False


@responses.activate
def test_failed_authentication_does_not_publish_singleton(mock_base_url, mock_api_url):
    responses.add(responses.GET, mock_api_url + "/auth/verify", json={}, status=200)
    with pytest.raises(ValueError, match="Invalid authentication response"):
        new("test-key", mock_base_url)
    assert not exists()
    responses.add(responses.GET, mock_api_url + "/auth/verify", json=PROFILE, status=200)
    c = new("test-key", mock_base_url)
    assert c.username == "test-user"
    reset()


@pytest.mark.skipif("fork" not in multiprocessing.get_all_start_methods(), reason="fork is unavailable")
@pytest.mark.parametrize("locked", [False, True])
@responses.activate
def test_fork_discards_inherited_client_and_lock(locked, mock_base_url, mock_api_url):
    from swanlab.sdk.internal.core_python import client as runtime

    responses.add(responses.GET, mock_api_url + "/auth/verify", json=PROFILE)
    owned = new("parent-key", mock_base_url)
    context = multiprocessing.get_context("fork")
    receive, send = context.Pipe(duplex=False)

    def child():
        send.send(runtime.exists())
        send.close()

    process = context.Process(target=child)
    try:
        if locked:
            runtime._client_lock.acquire()
        try:
            process.start()
        finally:
            if locked:
                runtime._client_lock.release()
        assert receive.poll(3), "child blocked on inherited client lock"
        assert receive.recv() is False
        process.join(3)
        assert process.exitcode == 0
        assert _get_client() is owned
    finally:
        if process.is_alive():
            process.terminate()
            process.join(3)
        receive.close()
        send.close()
        reset(owned)
