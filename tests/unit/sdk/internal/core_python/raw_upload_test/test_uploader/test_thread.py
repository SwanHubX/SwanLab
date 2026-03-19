"""
@author: caddiesnew
@file: test_thread.py
@time: 2026/3/19
@description: uploader 上传线程幂等性测试
"""

from queue import Queue
from unittest.mock import patch

import yaml
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal import adapter
from swanlab.core_python.uploader import UploadType
from swanlab.core_python.uploader.thread import RecordQueue, ThreadPool, UploadCollector


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

    writer.put((UploadType.SCALAR_METRIC, [make_scalar_record().SerializeToString()]))

    with patch("swanlab.core_python.uploader.upload.trace_metrics") as mock_trace:
        collector.callback(reader)
        collector.callback(reader)

    assert mock_trace.call_count == 1
    payload = mock_trace.call_args.args[1]

    assert mock_trace.call_args.args[0] == "/house/metrics"
    assert payload["type"] == "scalar"
    assert len(payload["metrics"]) == 1
    assert payload["metrics"][0]["key"] == "train/loss"
    assert payload["metrics"][0]["index"] == 3
    assert payload["metrics"][0]["type"] == "scalar"
    assert payload["metrics"][0]["scalar"] == {"number": 0.125}
    assert "timestamp" in payload["metrics"][0]


def test_threadpool_finish_is_idempotent_for_file_records(tmp_path):
    files_dir = tmp_path / adapter.dirname.files
    files_dir.mkdir()
    config = {"lr": 0.1, "epochs": 12}
    (files_dir / adapter.filename.config).write_text(yaml.safe_dump(config), encoding="utf-8")

    pool = ThreadPool(files_dir=files_dir)
    pool.put([make_config_record(), make_config_record()])

    with patch("swanlab.core_python.uploader.upload.trace_metrics") as mock_trace:
        pool.finish()
        pool.finish()

    assert mock_trace.call_count == 1
    assert mock_trace.call_args.args[0] == "/profile"
    assert mock_trace.call_args.args[1] == {"config": config}
    assert mock_trace.call_args.kwargs == {"method": "put", "per_request_len": -1}
