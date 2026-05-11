from pathlib import Path

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord, UpdateType
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord, ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext

# ── Default values matching CoreSettings ──

DEFAULT_BATCH_INTERVAL = 5.0
DEFAULT_BATCH_SIZE = 10_000


def _build_mock_ctx(**core_overrides):
    """Build a mock CoreContext with overridable core settings."""
    config = CoreConfig(
        run_id="test-run-id",
        run_dir=Path("/tmp/test-run"),
        section_rule=0,
        record_batch=core_overrides.get("batch_size", DEFAULT_BATCH_SIZE),
        record_interval=core_overrides.get("batch_interval", DEFAULT_BATCH_INTERVAL),
        save_split=100 * 1024 * 1024,
        save_size=50 * 1024 * 1024 * 1024,
        save_part=32 * 1024 * 1024,
        save_batch=100,
    )
    return CoreContext(config=config)


@pytest.fixture
def mock_ctx():
    """Mock CoreContext with default core settings. Tests read values from this instead of magic numbers."""
    return _build_mock_ctx()


@pytest.fixture
def make_ctx():
    """Factory fixture: create a mock CoreContext with overridable core settings."""

    def _make(**core_overrides):
        return _build_mock_ctx(**core_overrides)

    return _make


@pytest.fixture
def make_scalar_record():
    """工厂 fixture：创建一个 scalar 类型的 metric Record。"""

    def _make(step: int = 1) -> Record:
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        return Record(
            scalar=ScalarRecord(
                key="train/loss",
                step=step,
                timestamp=timestamp,
                type=ColumnType.COLUMN_TYPE_SCALAR,
                value=ScalarValue(number=0.125),
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
