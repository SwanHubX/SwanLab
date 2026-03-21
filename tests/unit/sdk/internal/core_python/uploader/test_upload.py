"""
@author: caddiesnew
@file: test_upload.py.py
@time: 2026/3/19 19:38
@description:
"""

from unittest.mock import MagicMock, call, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.upload import CoreTransportConfig, load_record, trace_records


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


def test_trace_records_upserts_every_record():
    transport = MagicMock()
    records = [make_scalar_record(step=1).SerializeToString(), make_scalar_record(step=2)]

    trace_records(records, per_request_len=-1, transport=transport)

    assert transport.upsert_record.call_count == 2
    first_record = transport.upsert_record.call_args_list[0].args[0]
    second_record = transport.upsert_record.call_args_list[1].args[0]
    assert first_record.metric.step == 1
    assert second_record.metric.step == 2
    transport.close.assert_not_called()


def test_trace_records_uploads_all_records_in_batches():
    transport = MagicMock()
    callback = MagicMock()
    records = [make_scalar_record(step=index) for index in range(2500)]

    with patch("swanlab.sdk.internal.core_python.uploader.upload.create_record_transport", return_value=transport) as factory:
        with patch("swanlab.sdk.internal.core_python.uploader.upload.time.sleep") as mock_sleep:
            trace_records(records, per_request_len=1000, upload_callback=callback)

    assert transport.upsert_record.call_count == 2500
    assert callback.call_args_list == [call(1000), call(1000), call(500)]
    assert mock_sleep.call_count == 3
    factory.assert_called_once_with()
    transport.close.assert_called_once_with()


def test_load_record_accepts_serialized_bytes():
    record = make_scalar_record(step=7)
    loaded = load_record(record.SerializeToString())
    assert loaded.metric.step == 7


def test_core_transport_config_reads_env(monkeypatch):
    monkeypatch.setenv("SWANLAB_CORE_HOST", "127.0.0.1")
    monkeypatch.setenv("SWANLAB_CORE_PORT", "9098")
    monkeypatch.setenv("SWANLAB_CORE_TIMEOUT", "5.5")

    config = CoreTransportConfig.from_env()

    assert config.enabled is True
    assert config.address == "127.0.0.1:9098"
    assert config.timeout == 5.5
