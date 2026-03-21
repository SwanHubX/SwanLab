"""
@author: caddiesnew
@file: test_upload.py.py
@time: 2026/3/19 19:38
@description: 
"""

from unittest.mock import MagicMock, call, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.batch import load_record, trace_records


@pytest.fixture
def mock_client():
    client = MagicMock()
    client._base_url = "https://mock.example.com/api"
    client._version = "1.2.3"
    client._before_request = MagicMock()
    client._session.cookies.get.return_value = "mock-sid"
    return client


def make_scalar_record(step: int = 1) -> Record:
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


def test_trace_records_upserts_every_record(mock_client):
    stub = MagicMock()
    records = [make_scalar_record(step=1).SerializeToString(), make_scalar_record(step=2)]

    with patch("swanlab.sdk.internal.core_python.uploader.batch.get_client", return_value=mock_client):
        with patch("swanlab.sdk.internal.core_python.uploader.batch._get_record_service", return_value=stub):
            trace_records(records, per_request_len=-1)

    assert stub.UpsertRecord.call_count == 2
    first_record = stub.UpsertRecord.call_args_list[0].args[0]
    second_record = stub.UpsertRecord.call_args_list[1].args[0]
    assert first_record.metric.step == 1
    assert second_record.metric.step == 2
    assert stub.UpsertRecord.call_args_list[0].kwargs["metadata"] == [
        ("cookie", "sid=mock-sid"),
        ("x-swanlab-sdk-version", "1.2.3"),
    ]
    assert stub.UpsertRecord.call_args_list[0].kwargs["wait_for_ready"] is True
    assert mock_client._before_request.call_count == 1


def test_trace_records_uploads_all_records_in_batches(mock_client):
    stub = MagicMock()
    callback = MagicMock()
    records = [make_scalar_record(step=index) for index in range(2500)]

    with patch("swanlab.sdk.internal.core_python.uploader.batch.get_client", return_value=mock_client):
        with patch("swanlab.sdk.internal.core_python.uploader.batch._get_record_service", return_value=stub):
            with patch("swanlab.sdk.internal.core_python.uploader.batch.time.sleep") as mock_sleep:
                trace_records(records, per_request_len=1000, upload_callback=callback)

    assert stub.UpsertRecord.call_count == 2500
    assert callback.call_args_list == [call(1000), call(1000), call(500)]
    assert mock_sleep.call_count == 3
    assert mock_client._before_request.call_count == 3


def test_load_record_accepts_serialized_bytes():
    record = make_scalar_record(step=7)
    loaded = load_record(record.SerializeToString())
    assert loaded.metric.step == 7
