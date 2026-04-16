import threading
from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch

# ─────────────────── group + route ───────────────────


def test_dispatch_groups_by_type(make_scalar_record, make_config_record):
    """混合类型 record 被正确分组分发。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer)

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(dispatch, "_upload_typed") as mock_upload:
        dispatch(metric_records + config_records)
        calls = mock_upload.call_args_list
        # metric 先出现所以先分发
        assert calls[0] == (("metric", metric_records),)
        assert calls[1] == (("config", config_records),)


def test_dispatch_calls_correct_handler(make_scalar_record):
    """各 _handle_{kind} 被调用且参数正确。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer)

    records = [make_scalar_record(step=1)]

    with patch.object(dispatch, "_upload_typed") as mock_upload:
        dispatch(records)
        mock_upload.assert_called_once_with("metric", records)


# ─────────────────── error rollback ───────────────────


def test_dispatch_error_rollback(make_scalar_record):
    """上传失败回滚到 buffer 头部（原地插入，buffer 对象引用不变）。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer)

    records = [make_scalar_record(step=1)]

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.create_record_sender") as mock_create:
        mock_sender = MagicMock()
        mock_sender.upload.side_effect = RuntimeError("upload failed")
        mock_create.return_value = mock_sender

        with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
            dispatch(records)

    # 回滚通过 self._buffer[:0] = records 原地插入，buffer 对象引用不变
    assert len(buffer) == 1
    assert buffer[0] is records[0]
    assert dispatch._buffer is buffer


# ─────────────────── unknown type ───────────────────


def test_dispatch_skips_unknown_type(make_scalar_record):
    """未知 kind 无 handler 时不报错，静默跳过。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer)

    # 手动构造一个不会匹配任何 _handle_{kind} 的场景
    # 通过 mock group_records_by_type 返回一个未知 key
    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"unknown_kind": [make_scalar_record(step=1)]})

        with patch.object(dispatch, "_upload_typed") as mock_upload:
            dispatch([make_scalar_record(step=1)])
            # unknown_kind 没有 _handle_unknown_kind，getattr 返回 None，不调用
            mock_upload.assert_not_called()
