"""
@author: caddiesnew
@file: test_upload.py.py
@time: 2026/3/19 19:38
@description:
"""

from unittest.mock import MagicMock, patch

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.sender import (
    NoopHttpRecordSender,
    create_http_record_sender,
    group_records_by_type,
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


def test_group_records_by_type_groups_records_by_record_type():
    buckets = group_records_by_type([make_scalar_record(step=1), make_config_record(), make_scalar_record(step=2)])

    assert set(buckets) == {"metric", "config"}
    assert [record.metric.step for record in buckets["metric"]] == [1, 2]
    assert len(buckets["config"]) == 1
    assert buckets["config"][0].config.update_type == UpdateType.UPDATE_TYPE_PATCH


def test_group_records_by_type_rejects_record_without_record_type():
    with pytest.raises(ValueError, match="record_type"):
        group_records_by_type([Record()])


def test_create_http_record_sender_returns_noop_sender():
    sender = create_http_record_sender()

    assert isinstance(sender, NoopHttpRecordSender)


def test_noop_http_record_sender_announces_once_per_instance():
    with patch("swanlab.sdk.internal.core_python.uploader.sender.console.debug") as mock_debug:
        sender_a = NoopHttpRecordSender()
        sender_b = NoopHttpRecordSender()

        sender_a.send(MagicMock())
        sender_a.send(MagicMock())
        sender_b.send(MagicMock())

    assert mock_debug.call_count == 2
