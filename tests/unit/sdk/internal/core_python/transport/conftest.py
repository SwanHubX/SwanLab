import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record


@pytest.fixture
def make_scalar_record():
    """工厂 fixture：创建一个 scalar 类型的 metric Record。"""

    def _make(step: int = 1) -> Record:
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

    return _make


@pytest.fixture
def make_config_record():
    """工厂 fixture：创建一个 config Record。"""

    def _make() -> Record:
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        return Record(config=ConfigRecord(update_type=UpdateType.UPDATE_TYPE_PATCH, timestamp=timestamp))

    return _make
