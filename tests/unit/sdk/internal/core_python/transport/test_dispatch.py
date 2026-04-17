import threading
from unittest.mock import MagicMock, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.dispatch import Dispatch

# ─────────────────── group + route ───────────────────


def test_dispatch_groups_by_type(make_scalar_record, make_config_record):
    """混合类型 record 被正确分组分发。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=set())

    metric_records = [make_scalar_record(step=1), make_scalar_record(step=2)]
    config_records = [make_config_record()]

    with patch.object(dispatch, "_handle_record_by_type", return_value=[]) as mock_handle:
        dispatch(metric_records + config_records)
        calls = mock_handle.call_args_list
        # metric 先出现所以先分发
        assert calls[0] == (("metric", metric_records),)
        assert calls[1] == (("config", config_records),)


def test_dispatch_calls_correct_handler(make_scalar_record):
    """各 _handle_{kind} 被调用且参数正确。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=set())

    records = [make_scalar_record(step=1)]

    with patch.object(dispatch, "_handle_record_by_type", return_value=[]) as mock_handle:
        dispatch(records)
        mock_handle.assert_called_once_with("metric", records)


# ─────────────────── error rollback ───────────────────


def test_dispatch_error_rollback(make_scalar_record):
    """上传失败回滚到 buffer 头部（原地插入，buffer 对象引用不变）。"""
    cond = threading.Condition()
    buffer = []
    buffer_num_index = set()

    mock_sender = MagicMock()
    mock_sender.upload.side_effect = RuntimeError("upload failed")
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=buffer_num_index, sender=mock_sender)

    records = [make_scalar_record(step=1)]
    records[0].num = 11

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        dispatch(records)

    # 回滚通过 self._buffer[:0] = records 原地插入，buffer 对象引用不变
    assert len(buffer) == 1
    assert buffer[0] is records[0]
    assert dispatch._buffer is buffer
    assert buffer_num_index == {11}


def test_dispatch_rollback_skips_records_already_buffered_by_num(make_scalar_record):
    """回滚时若 buffer 中已有同 num record，不再重复插入。"""
    cond = threading.Condition()
    existing = make_scalar_record(step=1)
    existing.num = 21
    buffer = [existing]
    buffer_num_index = {21}
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=buffer_num_index)

    duplicate = make_scalar_record(step=1)
    duplicate.num = 21
    other = make_scalar_record(step=2)
    other.num = 22

    with patch.object(dispatch, "_handle_record_by_type", return_value=[duplicate, other]):
        dispatch([duplicate, other])

    assert [record.num for record in buffer] == [22, 21]
    assert buffer_num_index == {21, 22}


# ─────────────────── unknown type ───────────────────


def test_dispatch_skips_unknown_type(make_scalar_record):
    """未知 kind 无 handler 时不报错，静默跳过。"""
    cond = threading.Condition()
    buffer = []
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=set())

    # 手动构造一个不会匹配任何 _handle_{kind} 的场景
    # 通过 mock group_records_by_type 返回一个未知 key
    with patch("swanlab.sdk.internal.core_python.transport.dispatch.group_records_by_type") as mock_group:
        from collections import OrderedDict

        mock_group.return_value = OrderedDict({"unknown_kind": [make_scalar_record(step=1)]})

        with patch.object(dispatch, "_upload_chunk", return_value=True) as mock_upload:
            dispatch([make_scalar_record(step=1)])
            # unknown_kind 在 _handle_record_by_type 内部被跳过，不会调用 _upload_chunk
            mock_upload.assert_not_called()


# ─────────────────── mixed type rollback ───────────────────


def test_dispatch_mixed_type_partial_failure_rollback(make_scalar_record, make_config_record):
    """混合类型中一种上传失败时，后续类型的 records 也被统一回滚到 buffer。"""
    cond = threading.Condition()
    buffer = []
    buffer_num_index = set()

    mock_sender = MagicMock()

    def upload_side_effect(record_type, records):
        if record_type == "metric":
            raise RuntimeError("upload failed")

    mock_sender.upload.side_effect = upload_side_effect
    dispatch = Dispatch(cond=cond, buffer=buffer, buffer_num_index=buffer_num_index, sender=mock_sender)

    metric_records = [make_scalar_record(step=1)]
    metric_records[0].num = 31
    config_records = [make_config_record()]
    config_records[0].num = 32

    with patch("swanlab.sdk.internal.core_python.transport.dispatch.console"):
        dispatch(metric_records + config_records)

    # metric 失败 + config 未处理，两者都应回滚到 buffer
    assert len(buffer) == 2
    assert buffer[0] is metric_records[0]
    assert buffer[1] is config_records[0]
    assert buffer_num_index == {31, 32}


def test_dispatch_record_types_follow_proto_descriptor():
    """分发类型集合应直接来自 Record.record_type oneof。"""
    assert Dispatch._RECORD_TYPES == frozenset(f.name for f in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)
