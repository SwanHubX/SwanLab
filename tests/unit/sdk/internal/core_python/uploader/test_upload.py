"""
@author: caddiesnew
@file: test_upload.py.py
@time: 2026/3/19 19:38
@description:
"""

from typing import Sequence, cast
from unittest.mock import MagicMock, call, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.sender import (
    CoreTransportConfig,
    HttpRecordTransport,
    create_record_transport,
    trace_records,
)


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


def make_config_record() -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(config=ConfigRecord(update_type=UpdateType.UPDATE_TYPE_PATCH, timestamp=timestamp))


def test_trace_records_groups_records_by_record_type():
    transport = MagicMock()
    metric_1 = make_scalar_record(step=1)
    config = make_config_record()
    metric_2 = make_scalar_record(step=2)
    records = [metric_1, config, metric_2]

    trace_records(records, per_request_len=-1, transport=transport)

    assert transport.upload_record_group.call_args_list == [
        call("metric", [metric_1, metric_2]),
        call("config", [config]),
    ]
    transport.close.assert_not_called()


def test_trace_records_uploads_all_records_in_batches():
    transport = MagicMock()
    callback = MagicMock()
    records = [make_scalar_record(step=index) for index in range(2500)]

    with patch(
        "swanlab.sdk.internal.core_python.uploader.sender.create_record_transport", return_value=transport
    ) as factory:
        with patch("swanlab.sdk.internal.core_python.uploader.sender.time.sleep") as mock_sleep:
            trace_records(records, per_request_len=1000, upload_callback=callback)

    assert transport.upload_record_group.call_args_list == [
        call("metric", records[:1000]),
        call("metric", records[1000:2000]),
        call("metric", records[2000:]),
    ]
    assert callback.call_args_list == [call(1000), call(1000), call(500)]
    assert mock_sleep.call_count == 3
    factory.assert_called_once_with()
    transport.close.assert_called_once_with()


def test_trace_records_rejects_serialized_bytes():
    transport = MagicMock()
    bad_records = cast(Sequence[Record], [make_scalar_record(step=7).SerializeToString()])

    try:
        trace_records(bad_records, transport=transport)
    except TypeError as exc:
        assert "Record" in str(exc)
    else:
        raise AssertionError("trace_records should reject serialized bytes")


def test_create_record_transport_returns_http_transport_when_enabled():
    transport = create_record_transport(CoreTransportConfig(enabled=True, address=None, timeout=1.0))

    assert isinstance(transport, HttpRecordTransport)


def test_http_record_transport_dispatches_record_group_to_typed_handler():
    transport = HttpRecordTransport(CoreTransportConfig(enabled=True, address=None, timeout=1.0))
    records = [make_scalar_record(step=1)]

    with patch.object(transport, "_upload_metric_group") as mock_handler:
        transport.upload_record_group("metric", records)

    mock_handler.assert_called_once_with(records)


def test_core_transport_config_reads_env(monkeypatch):
    monkeypatch.setenv("SWANLAB_CORE_HOST", "127.0.0.1")
    monkeypatch.setenv("SWANLAB_CORE_PORT", "9098")
    monkeypatch.setenv("SWANLAB_CORE_TIMEOUT", "5.5")

    config = CoreTransportConfig.from_env()

    assert config.enabled is True
    assert config.address == "127.0.0.1:9098"
    assert config.timeout == 5.5
