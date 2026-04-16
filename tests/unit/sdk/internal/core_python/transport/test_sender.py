from typing import Sequence, cast

from swanlab.sdk.internal.core_python.transport.sender import (
    HttpRecordTransport,
    create_record_transport,
)


def test_http_record_transport_upload_record_group_skips_empty():
    """验证上传空列表时不抛异常，直接跳过。"""
    transport = HttpRecordTransport()
    transport.upload_record_group("metric", [])  # should not raise


def test_http_record_transport_close():
    """验证 close() 不抛异常（当前为空实现）。"""
    transport = HttpRecordTransport()
    transport.close()  # should not raise


# ─────────────────── create_record_transport ───────────────────


def test_create_record_transport_returns_http_transport():
    """验证工厂函数返回 HttpRecordTransport 实例。"""
    transport = create_record_transport()

    assert isinstance(transport, HttpRecordTransport)
