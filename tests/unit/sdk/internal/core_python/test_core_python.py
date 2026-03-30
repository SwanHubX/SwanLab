from unittest.mock import MagicMock, patch

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.metric.data.v1.scalar.scalar_pb2 import ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunConfig, RunContext
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.settings import Settings


def make_scalar_record() -> Record:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return Record(
        metric=DataRecord(
            key="train/loss",
            step=1,
            timestamp=timestamp,
            type=ColumnType.COLUMN_TYPE_FLOAT,
            scalar=ScalarValue(number=0.125),
        )
    )


def make_ctx(tmp_path) -> RunContext:
    settings = Settings.model_validate({"mode": "cloud", "run": {"id": "test-run-id"}})
    return RunContext(RunConfig(run_dir=tmp_path, settings=settings))


def test_core_python_cloud_mode_forwards_records_to_uploader(tmp_path):
    ctx = make_ctx(tmp_path)

    with patch("swanlab.sdk.internal.core_python.Uploader") as mock_pool_cls:
        uploader = MagicMock()
        mock_pool_cls.return_value = uploader

        core = CorePython(ctx)
        core.startup(cloud=True, persistence=False)

        record = make_scalar_record()
        core.handle_records([record])
        uploader.put.assert_called_once_with([record])

        core.shutdown()
        uploader.finish.assert_called_once_with()
