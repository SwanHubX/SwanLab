"""
@author: caddiesnew
@file: test_thread.py
@time: 2026/3/19
@description: uploader 上传线程幂等性测试
"""

import time
from queue import Queue
from unittest.mock import patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.batch import load_record
from swanlab.sdk.internal.core_python.uploader.thread import RecordQueue, ThreadPool, UploadCollector


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
    queue = Queue()
    writer = RecordQueue(queue=queue, readable=False, writable=True)
    reader = RecordQueue(queue=queue, readable=True, writable=False)
    collector = UploadCollector()

    writer.put(make_scalar_record().SerializeToString())

    with patch("swanlab.sdk.internal.core_python.uploader.thread.log_collector.upload_records") as mock_upload:
        collector.callback(reader)
        collector.callback(reader)

    assert mock_upload.call_count == 1
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 1
    record = load_record(uploaded[0])
    assert record.metric.key == "train/loss"
    assert record.metric.step == 3
    assert record.metric.scalar.number == 0.125


def test_threadpool_finish_is_idempotent_for_config_records():
    pool = ThreadPool()
    pool.put([make_config_record(), make_config_record()])

    with patch("swanlab.sdk.internal.core_python.uploader.thread.log_collector.upload_records") as mock_upload:
        pool.finish()
        pool.finish()

    assert mock_upload.call_count == 1
    uploaded = mock_upload.call_args.args[0]
    assert len(uploaded) == 2
    assert all(load_record(record).WhichOneof("record_type") == "config" for record in uploaded)


def test_threadpool_starts_upload_thread_automatically():
    with patch.object(ThreadPool, "SLEEP_TIME", 0.01):
        pool = ThreadPool(upload_interval=0.01)

        try:
            with patch("swanlab.sdk.internal.core_python.uploader.thread.log_collector.upload_records") as mock_upload:
                pool.put([make_scalar_record()])

                deadline = time.time() + 1.0
                while mock_upload.call_count == 0 and time.time() < deadline:
                    time.sleep(0.02)

                assert mock_upload.call_count == 1
        finally:
            pool.finish()
