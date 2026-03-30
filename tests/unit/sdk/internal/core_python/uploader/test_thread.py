"""
@author: caddiesnew
@file: test_thread.py
@time: 2026/3/19
@description: uploader 上传线程幂等性测试
"""

import threading
import time
from unittest.mock import MagicMock, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.bus.events import MetricsUploadEvent
from swanlab.sdk.internal.core_python.uploader.helper import RecordQueue
from swanlab.sdk.internal.core_python.uploader.uploader import HttpBatchUploader


def make_scalar_record(step: int = 3) -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        metric=DataRecord(
            key="train/loss",
            step=step,
            timestamp=timestamp,
            type=ColumnType.COLUMN_TYPE_FLOAT,
            scalar=ScalarValue(number=0.125),
        )
    )


def make_config_record() -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(config=ConfigRecord(update_type=UpdateType.UPDATE_TYPE_PATCH, timestamp=timestamp))


def test_record_queue_put_and_get_roundtrip():
    queue = RecordQueue()
    record = make_scalar_record()

    queue.put(record)

    loaded = queue.get_nowait()
    assert loaded.metric.key == "train/loss"
    assert queue.empty()


def test_http_batch_uploader_flush_groups_records_by_type():
    sender = MagicMock()
    uploader = HttpBatchUploader(sender=sender, auto_start=False)
    uploader.enqueue([make_config_record(), make_scalar_record(), make_scalar_record()])

    uploader.flush()

    assert sender.send.call_count == 1
    upload_event = sender.send.call_args.args[0]
    assert isinstance(upload_event, MetricsUploadEvent)
    buckets = upload_event.buckets
    assert set(buckets) == {"config", "metric"}
    assert len(buckets["config"]) == 1
    assert len(buckets["metric"]) == 2
    assert uploader._pending == []
    assert uploader._queue.empty()


def test_http_batch_uploader_flush_keeps_pending_records_when_send_fails():
    sender = MagicMock()
    sender.send.side_effect = RuntimeError("send boom")
    uploader = HttpBatchUploader(sender=sender, auto_start=False)
    uploader.enqueue([make_config_record(), make_scalar_record()])

    with pytest.raises(RuntimeError, match="send boom"):
        uploader.flush()

    assert len(uploader._pending) == 2
    assert uploader._queue.empty()


def test_http_batch_uploader_drops_oldest_records_when_pending_limit_is_exceeded():
    sender = MagicMock()
    sender.send.side_effect = RuntimeError("send boom")
    uploader = HttpBatchUploader(sender=sender, auto_start=False, max_pending=3)
    uploader.enqueue([make_scalar_record(step=1), make_scalar_record(step=2)])

    with pytest.raises(RuntimeError, match="send boom"):
        uploader.flush()

    uploader.enqueue([make_scalar_record(step=3), make_scalar_record(step=4), make_scalar_record(step=5)])

    with patch("swanlab.sdk.internal.core_python.uploader.uploader.console.warning") as mock_warning:
        with pytest.raises(RuntimeError, match="send boom"):
            uploader.flush()

    assert [record.metric.step for record in uploader._pending] == [3, 4, 5]
    mock_warning.assert_called_once()


def test_http_batch_uploader_close_flushes_pending_records_and_closes_sender():
    sender = MagicMock()
    uploader = HttpBatchUploader(sender=sender, auto_start=False)
    uploader.enqueue([make_scalar_record()])

    uploader.close()

    assert sender.send.call_count == 1
    sender.close.assert_called_once_with()


def test_http_batch_uploader_close_closes_sender_when_final_flush_raises():
    sender = MagicMock()
    sender.send.side_effect = RuntimeError("flush boom")
    uploader = HttpBatchUploader(sender=sender, auto_start=False)
    uploader.enqueue([make_scalar_record()])

    with pytest.raises(RuntimeError, match="flush boom"):
        uploader.close()

    sender.close.assert_called_once_with()


def test_http_batch_uploader_starts_timer_and_flushes_automatically():
    sender = MagicMock()
    uploader = HttpBatchUploader(sender=sender, upload_interval=0.01)

    try:
        uploader.enqueue([make_scalar_record()])

        deadline = time.time() + 1.0
        while sender.send.call_count == 0 and time.time() < deadline:
            time.sleep(0.02)

        assert sender.send.call_count == 1
    finally:
        uploader.close()


def test_http_batch_uploader_close_cancels_and_joins_timer():
    sender = MagicMock()

    with patch("swanlab.sdk.internal.core_python.uploader.uploader.Timer") as mock_timer_cls:
        timer = MagicMock()
        mock_timer_cls.return_value = timer

        uploader = HttpBatchUploader(sender=sender)
        uploader.close()

    timer.start.assert_called_once_with()
    timer.cancel.assert_called_once_with()
    timer.join.assert_called_once_with(timeout=10)


def test_http_batch_uploader_rejects_enqueue_after_close():
    uploader = HttpBatchUploader(sender=MagicMock(), auto_start=False)
    uploader.close()

    with pytest.raises(RuntimeError, match="closed"):
        uploader.enqueue([make_scalar_record()])


def test_http_batch_uploader_close_waits_for_inflight_enqueue_before_final_flush():
    sender = MagicMock()
    uploader = HttpBatchUploader(sender=sender, auto_start=False)
    put_started = threading.Event()
    allow_put = threading.Event()
    original_put_all = uploader._queue.put_all

    def blocking_put_all(records):
        put_started.set()
        allow_put.wait(timeout=1)
        original_put_all(records)

    uploader._queue.put_all = blocking_put_all

    enqueue_thread = threading.Thread(target=uploader.enqueue, args=([make_scalar_record()],))
    enqueue_thread.start()
    assert put_started.wait(timeout=1)

    close_thread = threading.Thread(target=uploader.close)
    close_thread.start()
    time.sleep(0.05)
    assert close_thread.is_alive()

    allow_put.set()
    enqueue_thread.join(timeout=1)
    close_thread.join(timeout=1)

    assert sender.send.call_count == 1
    sender.close.assert_called_once_with()
