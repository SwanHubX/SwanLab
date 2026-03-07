"""
@author: cunyue
@file: test_client_session.py
@time: 2026/3/7 21:32
@description: 测试 SwanLab 运行时客户端会话辅助函数
"""

import pytest
import responses

from swanlab.sdk.internal.core_python.client.session import TimeoutHTTPAdapter, create, request_retries_ctx
from swanlab.sdk.pkg.exceptions import ApiError


@pytest.fixture
def session():
    """创建一个默认的 Session 实例用于测试"""
    return create(timeout=10, default_retry=0)  # 为了测试报错包装，关掉默认重试


def test_create_session_properties():
    """测试会话创建时的默认属性挂载"""
    sess = create(timeout=30)
    # 验证是否挂载了 Version Header
    assert "X-SwanLab-SDK-Version" in sess.headers

    sess = create(timeout=30, default_retry=5)
    # 验证是否挂载了 Version Header
    assert "X-SwanLab-SDK-Version" in sess.headers

    # 获取 adapter
    adapter = sess.get_adapter("http://")

    # 加入类型断言：告诉类型检查器这具体是一个 TimeoutHTTPAdapter
    assert isinstance(adapter, TimeoutHTTPAdapter)

    # 现在的 adapter 会被 IDE 识别为 TimeoutHTTPAdapter，不再报错
    assert adapter.timeout == 30
    assert adapter.max_retries.total == 5


@responses.activate()
def test_session_send_success(session):
    """测试 2xx 成功响应时，直接返回 Response 对象"""
    responses.add(responses.GET, "http://api.test/data", json={"ok": True}, status=200)

    resp = session.get("http://api.test/data")
    assert resp.status_code == 200
    assert resp.json() == {"ok": True}


@responses.activate()
def test_session_send_api_error_custom_parsed(session):
    """测试 4xx/5xx 时，能否正确提取 JSON 中的错误信息并抛出 ApiError"""
    responses.add(
        responses.POST,
        "http://api.test/submit",
        json={"code": "1001", "message": "API Token Invalid"},
        status=401,
        headers={"traceid": "trace-abc-123"},
    )

    with pytest.raises(ApiError) as exc_info:
        session.post("http://api.test/submit")

    error = exc_info.value
    assert error.code == "1001"
    assert error.message == "API Token Invalid"
    assert error.trace_id == "trace-abc-123"
    assert error.response.status_code == 401


@responses.activate()
def test_session_send_api_error_fallback(session):
    """测试 502 等非 JSON 报错时，能否回退为通用报错信息抛出"""
    responses.add(responses.GET, "http://api.test/bad", body="<html>502 Bad Gateway</html>", status=502)

    with pytest.raises(ApiError) as exc_info:
        session.get("http://api.test/bad")

    error = exc_info.value
    # 因为没解析出 code，所以默认会设置成未知或回退值
    assert error.code == "unknown code"
    assert error.message == "unknown error"
    assert error.trace_id == "unknown"  # 没有 traceid header


@responses.activate()
def test_request_retries_context_cleanup(session):
    """测试无论请求成功与否，重试次数的 ContextVar 都会被清理"""
    responses.add(responses.GET, "http://api.test/fail", status=500)

    # 确认初始状态为空
    assert request_retries_ctx.get() is None

    # 带 retries 参数发起请求，预期会失败
    with pytest.raises(ApiError):
        session.get("http://api.test/fail", retries=3)

    # 无论中间发生了什么（包括抛出异常），finally 块必须确保清理了上下文
    assert request_retries_ctx.get() is None
