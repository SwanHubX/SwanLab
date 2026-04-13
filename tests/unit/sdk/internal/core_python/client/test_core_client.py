"""
@author: cunyue
@file: test_core_client.py
@time: 2026/3/7 20:45
@description: 测试 SwanLab 运行时客户端
"""

from datetime import datetime, timedelta, timezone
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.core_python.client import (
    _get_client,
    exists,
    new,
    reset,
)
from swanlab.sdk.internal.core_python.client import get as global_get


@pytest.fixture()
def mock_base_url():
    """用户传入的基础 URL"""
    return "http://mock.example.com"


@pytest.fixture()
def mock_api_url(mock_base_url):
    """Client 内部强制拼接后的真实 API URL"""
    return f"{mock_base_url}/api"


@pytest.fixture()
def mock_login():
    patcher = patch("swanlab.sdk.internal.pkg.client.login_by_api_key")
    mock_func = patcher.start()
    mock_func.return_value = {"sid": "mock-token-123", "expiredAt": "2999-01-01T00:00:00.000Z"}
    yield mock_func
    patcher.stop()


# -------------------------------------------------------------------
# Test Cases (测试用例)
# -------------------------------------------------------------------


def test_global_proxy_functions(mock_login, mock_base_url, mock_api_url):
    """测试 new, exists, reset 全局状态的生命周期"""
    _ = mock_login

    # 清理初始状态以防干扰
    if exists():
        reset()

    # 1. 未初始化时调用，应该报错
    with pytest.raises(RuntimeError, match="not initialized"):
        _get_client()

    # 2. 正确初始化
    global_client = new("global-key", mock_base_url)
    assert exists() is True
    assert global_client._api_key == "global-key"

    # 3. 重复初始化应该报错
    with pytest.raises(RuntimeError, match="already exists"):
        new("another-key", mock_base_url)

    # 4. 测试全局快捷方法 (如 swanlab.client.get)
    global_client._expired_at = datetime.now(timezone.utc) + timedelta(days=30)
    global_client._session.request = MagicMock()
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"ok": True}
    global_client._session.request.return_value = mock_resp

    resp = global_get("/global-test")

    # 注意这里断言的应该是拼接后的 mock_api_url
    global_client._session.request.assert_called_with("GET", mock_api_url + "/global-test", params=None, retries=None)
    assert resp.data == {"ok": True}

    # 5. 销毁并验证
    reset()
    assert exists() is False
    with pytest.raises(RuntimeError, match="not initialized"):
        reset()  # 再次销毁应报错
