"""
@author: caddiesnew
@file: test_thread.py
@time: 2026/3/19
@description: uploader 上传线程幂等性测试
"""

import threading
import time
from unittest.mock import MagicMock, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader import Collector, ThreadPool


def make_scalar_record() -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        metric=DataRecord(
            key="train/loss",
            step=3,
            timestamp=timestamp,
            type=ColumnType.COLUMN_TYPE_FLOAT,
            scalar=ScalarValue(number=0.125),
        )
    )


def make_config_record() -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(config=ConfigRecord(update_type=UpdateType.UPDATE_TYPE_PATCH, timestamp=timestamp))


def test_upload_collector_callback_is_idempotent_for_scalar_records():
    collector = Collector()
    record = make_scalar_record()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        collector.callback([record])
        collector.callback([])

    assert mock_upload.call_count == 1
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 1
    record = uploaded[0]
    assert record.metric.key == "train/loss"
    assert record.metric.step == 3
    assert record.metric.scalar.number == 0.125


def test_upload_collector_callback_waits_for_lock_release():
    collector = Collector()
    record = make_scalar_record()
    done = threading.Event()

    collector._lock.acquire()

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        worker = threading.Thread(target=lambda: (collector.callback([record]), done.set()))
        worker.start()

        time.sleep(0.05)
        assert done.is_set() is False
        mock_upload.assert_not_called()

        collector._lock.release()
        worker.join(timeout=1)

    assert done.is_set() is True


def test_threadpool_drain_records_returns_pending_records_once():
    pool = ThreadPool(auto_start=False)
    records = [make_config_record(), make_config_record()]
    pool.put(records)

    assert pool._drain_records() == records
    assert pool._drain_records() == []


def test_threadpool_finish_is_idempotent_for_config_records():
    pool = ThreadPool(auto_start=False)
    pool.put([make_config_record(), make_config_record()])

    with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
        pool.finish()
        pool.finish()

    assert mock_upload.call_count == 1
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 2
    assert all(record.WhichOneof("record_type") == "config" for record in uploaded)


def test_threadpool_starts_upload_thread_automatically():
    with patch.object(ThreadPool, "UPLOAD_INTERVAL", 0.01):
        pool = ThreadPool(upload_interval=0.01)

        try:
            with patch("swanlab.sdk.internal.core_python.uploader.collector.upload_records") as mock_upload:
                pool.put([make_scalar_record()])

                deadline = time.time() + 1.0
                while mock_upload.call_count == 0 and time.time() < deadline:
                    time.sleep(0.02)

                assert mock_upload.call_count == 1
        finally:
            pool.finish()


def test_threadpool_start_reuses_pkg_timer_scheduler():
    timer = MagicMock()

    with patch.object(ThreadPool, "UPLOAD_INTERVAL", 0.25), patch(
        "swanlab.sdk.internal.core_python.uploader.thread.Timer"
    ) as mock_timer_cls:
        mock_timer_cls.return_value = timer

        pool = ThreadPool(auto_start=False)
        pool.start()

        assert mock_timer_cls.call_count == 1
        _, kwargs = mock_timer_cls.call_args
        assert kwargs["interval"] == 0.25
        assert kwargs["immediate"] is True
        assert kwargs["name"] == ThreadPool.UPLOAD_THREAD_NAME
        timer.start.assert_called_once_with()


def test_threadpool_finish_cancels_and_joins_timer():
    timer = MagicMock()

    with patch("swanlab.sdk.internal.core_python.uploader.thread.Timer") as mock_timer_cls:
        mock_timer_cls.return_value = timer

        pool = ThreadPool(auto_start=False)
        pool.start()
        pool.finish()

    timer.cancel.assert_called_once_with()
    timer.join.assert_called_once_with(timeout=10)
