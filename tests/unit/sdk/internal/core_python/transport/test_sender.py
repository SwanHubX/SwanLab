from unittest.mock import call, patch

from swanlab.sdk.internal.core_python.transport.sender import (
    HttpRecordSender,
    create_record_sender,
)

# ─────────────────── HttpRecordSender.upload ───────────────────


def test_upload_skips_empty():
    """验证上传空列表时不调用任何 upload_{kind} 方法。"""
    sender = HttpRecordSender()
    with patch.object(sender, "upload_metric") as mock_fn:
        sender.upload("metric", [])
        mock_fn.assert_not_called()


def test_upload_routes_to_correct_method(make_scalar_record):
    """验证 upload() 通过 getattr 路由到对应 upload_{kind} 方法。"""
    sender = HttpRecordSender()
    records = [make_scalar_record(step=1)]

    with patch.object(sender, "upload_metric") as mock_fn:
        sender.upload("metric", records)
        mock_fn.assert_called_once_with(records)


def test_upload_unknown_type_skipped(make_scalar_record):
    """验证 upload() 对未知 record_type 静默跳过（无对应方法时不报错）。"""
    sender = HttpRecordSender()
    records = [make_scalar_record(step=1)]
    # "unknown" 不存在 upload_unknown 方法，应静默跳过
    sender.upload("unknown", records)  # should not raise


# ─────────────────── upload_{kind} methods ───────────────────


def test_upload_kind_methods_are_callable():
    """验证所有 9 种 upload_{kind} 方法可调用且不抛异常。"""
    sender = HttpRecordSender()
    for kind in ("run", "finish", "column", "metric", "config", "console", "metadata", "requirements", "conda"):
        fn = getattr(sender, f"upload_{kind}")
        fn([])  # should not raise


def test_upload_metric_called_with_records(make_scalar_record):
    """验证 upload_metric 被传入正确的 records。"""
    sender = HttpRecordSender()
    records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    with patch("swanlab.sdk.internal.pkg.console.debug") as mock_debug:
        sender.upload_metric(records)
        mock_debug.assert_called_once()


# ─────────────────── close / factory ───────────────────


def test_http_record_sender_close():
    """验证 close() 不抛异常（当前为空实现）。"""
    sender = HttpRecordSender()
    sender.close()  # should not raise


def test_create_record_sender_returns_http_sender():
    """验证工厂函数返回 HttpRecordSender 实例。"""
    sender = create_record_sender()
    assert isinstance(sender, HttpRecordSender)
