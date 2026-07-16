import asyncio
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord, ScalarValue
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, ResumeMode, RunState, StartRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.core_python.utils import PrepareExperimentStartResult


def make_timestamp() -> Timestamp:
    ts = Timestamp()
    ts.GetCurrentTime()
    return ts


def make_core_config(run_dir: Path, run_id: str = "sync-run-id") -> CoreConfig:
    return CoreConfig(
        run_id=run_id,
        run_dir=run_dir,
        section_rule=0,
        record_batch=10000,
        record_interval=5.0,
        save_split=100 * 1024 * 1024,
        save_size=50 * 1024 * 1024 * 1024,
        save_part=32 * 1024 * 1024,
        save_batch=100,
    )


def make_core_settings(run_dir: Path, run_id: str = "sync-run-id") -> CoreSettings:
    config = make_core_config(run_dir, run_id)
    return CoreSettings(
        run_id=config.run_id,
        run_dir=str(config.run_dir),
        section_rule=config.section_rule,
        record_batch=config.record_batch,
        record_interval=config.record_interval,
        save_split=config.save_split,
        save_size=config.save_size,
        save_part=config.save_part,
        save_batch=config.save_batch,
    )


def make_start_record(**kwargs) -> StartRecord:
    values = {
        "id": "original-run-id",
        "workspace": "original-workspace",
        "project": "original-project",
        "name": "original-name",
        "color": "#123456",
        "started_at": make_timestamp(),
    }
    values.update(kwargs)
    return StartRecord(**values)


def make_prepare_result(new_experiment: bool = True) -> PrepareExperimentStartResult:
    return PrepareExperimentStartResult(
        username="alice",
        project="demo",
        project_data={
            "cuid": "project-id",
            "name": "demo",
            "group": {"username": "alice"},
            "path": "/alice/demo",
            "visibility": "PRIVATE",
            "_count": {"experiments": 1, "contributors": 0, "collaborators": 0, "clones": 0},
        },
        experiment={"cuid": "experiment-id", "slug": "sync-run-id", "name": "synced-run"},
        new_experiment=new_experiment,
        name="synced-run",
        color="#abcdef",
    )


def write_run_file(run_dir: Path, *records: Record) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / "run-sync-run-id.swanlab"
    writer = DataStoreWriter()
    writer.open(path)
    for record in records:
        writer.write(record.SerializeToString())
    writer.close()
    return path


def make_start_request(run_dir: Path, **kwargs) -> DeliverSyncStartRequest:
    return DeliverSyncStartRequest(core_settings=make_core_settings(run_dir), **kwargs)


def record_kind(record: Record) -> str:
    for kind in [
        "start",
        "finish",
        "column",
        "scalar",
        "media",
        "log",
        "save",
    ]:
        if record.HasField(kind):
            return kind
    return "unknown"


class FakeTransport:
    def __init__(self, ctx: CoreContext, tracker=None):
        self.ctx = ctx
        self.tracker = tracker
        self.records: list[Record] = []
        self.started = False
        self.finished = False
        self.joined = False

    def put(self, records: list[Record]) -> None:
        self.records.extend(records)

    def start(self) -> None:
        self.started = True

    def request_finish(self) -> None:
        self.finished = True
        if self.tracker is not None:
            self.tracker.set_state(CoreState.CORE_STATE_FINISHED)

    def join(self, timeout=None) -> bool:
        self.joined = True
        if self.tracker is not None:
            self.tracker.set_state(CoreState.CORE_STATE_FINISHED)
        return True

    def finish(self, timeout=None) -> bool:
        self.request_finish()
        return self.join(timeout=timeout)


class FakeExecutor:
    def __init__(self):
        self.started = False
        self.waited = False
        self.closed = False
        self.coro = None

    def start(self, coro):
        self.started = True
        self.coro = coro
        if hasattr(coro, "close"):
            coro.close()

    def wait(self):
        self.waited = True

    def close(self):
        self.closed = True


class FakeReader:
    def __init__(self, records: list[Record]):
        self._records = records
        self.closed = False

    def __iter__(self):
        return iter([record.SerializeToString() for record in self._records])

    def close(self):
        self.closed = True


class FakeRawReader:
    def __init__(self, payloads: list[bytes]):
        self._payloads = payloads
        self.closed = False

    def __iter__(self):
        return iter(self._payloads)

    def close(self):
        self.closed = True


class FakeMetric:
    def __init__(self):
        self.updated = []

    def ensure_type_match(self, _: ColumnType) -> None:
        pass

    def try_accept_step(self, _: int) -> bool:
        return True

    def update(self, record) -> None:
        self.updated.append(record)


class FakeMetrics:
    def __init__(self):
        self.metrics: dict[str, FakeMetric] = {}

    def get(self, key: str) -> Optional[FakeMetric]:
        return self.metrics.get(key)

    def define_scalar(self, **kwargs) -> FakeMetric:
        metric = FakeMetric()
        self.metrics[kwargs["key"]] = metric
        return metric

    def define_media(self, **kwargs) -> FakeMetric:
        metric = FakeMetric()
        self.metrics[kwargs["key"]] = metric
        return metric


def test_sync_context_run_file_finds_single_file(tmp_path: Path):
    run_file = tmp_path / "run-abc.swanlab"
    run_file.touch()
    ctx = CoreContext(config=make_core_config(tmp_path, run_id=""), mode="sync")

    assert ctx.run_file == run_file


def test_sync_context_run_file_missing_raises(tmp_path: Path):
    ctx = CoreContext(config=make_core_config(tmp_path, run_id=""), mode="sync")

    with pytest.raises(FileNotFoundError, match=r"No run-\*\.swanlab"):
        _ = ctx.run_file


def test_sync_context_run_file_multiple_raises(tmp_path: Path):
    (tmp_path / "run-a.swanlab").touch()
    (tmp_path / "run-b.swanlab").touch()
    ctx = CoreContext(config=make_core_config(tmp_path, run_id=""), mode="sync")

    with pytest.raises(RuntimeError, match=r"Multiple run-\*\.swanlab"):
        _ = ctx.run_file


def test_deliver_sync_start_reads_start_record(tmp_path: Path):
    start_record = make_start_record(id="start-id")
    write_run_file(tmp_path, Record(start=start_record))
    core = CoreSyncPython()

    resp = core.deliver_sync_start(make_start_request(tmp_path))

    assert resp.success is True
    assert core._start_record is not None
    assert core._start_record.id == "start-id"


def test_deliver_sync_start_rejects_missing_start_record(tmp_path: Path):
    write_run_file(
        tmp_path,
        Record(
            save=SaveRecord(
                name="config",
                source_path=(tmp_path / "config.yaml").as_posix(),
                type=SaveType.SAVE_TYPE_CONFIG,
            )
        ),
    )
    core = CoreSyncPython()

    resp = core.deliver_sync_start(make_start_request(tmp_path))

    assert resp.success is False
    assert "missing the required start record" in resp.message


def test_deliver_sync_flush_prepares_experiment_with_overrides(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._start_record = make_start_record()
    core._start_request = make_start_request(
        tmp_path,
        workspace="override-workspace",
        project="override-project",
        id="override-run-id",
    )
    prepare = MagicMock(return_value=make_prepare_result(new_experiment=True))
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.prepare_experiment_start", prepare)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.RunMetrics.new", MagicMock(return_value=(FakeMetrics(), -1, -1, -1))
    )
    transports: list[FakeTransport] = []

    def make_transport(ctx: CoreContext, tracker=None) -> FakeTransport:
        transport = FakeTransport(ctx, tracker=tracker)
        transports.append(transport)
        return transport

    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.Transport", make_transport)
    fake_executor = FakeExecutor()
    core._read_executor = fake_executor  # type: ignore[assignment]

    resp = core.deliver_sync_flush()

    assert resp.success is True
    assert resp.path == "/alice/demo/sync-run-id"
    prepared_record = prepare.call_args.args[0]
    assert prepared_record.workspace == "override-workspace"
    assert prepared_record.project == "override-project"
    assert prepared_record.id == "override-run-id"
    assert prepared_record.resume == ResumeMode.RESUME_MODE_ALLOW
    assert prepared_record.color == ""
    assert core._ctx.username == "alice"
    assert core._ctx.project == "demo"
    assert core._ctx.project_id == "project-id"
    assert core._ctx.experiment_id == "experiment-id"
    assert transports[0].tracker is core._tracker
    assert transports[0].started is True
    assert fake_executor.started is True
    assert core.get_operation_stats().stats.state == CoreState.CORE_STATE_RUNNING


def test_deliver_sync_flush_fetches_summary_for_existing_experiment(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._start_record = make_start_record()
    core._start_request = make_start_request(tmp_path)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start",
        MagicMock(return_value=make_prepare_result(new_experiment=False)),
    )
    get_summary = MagicMock(return_value={"log": None, "media": None, "scalar": None})
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.get_experiment_summary", get_summary)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.RunMetrics.new", MagicMock(return_value=(FakeMetrics(), -1, -1, -1))
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.Transport", lambda ctx, tracker=None: FakeTransport(ctx, tracker)
    )
    core._read_executor = FakeExecutor()  # type: ignore[assignment]

    resp = core.deliver_sync_flush()

    assert resp.success is True
    get_summary.assert_called_once_with("project-id", "experiment-id")


def test_deliver_sync_flush_returns_failure_when_prepare_fails(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._start_record = make_start_record()
    core._start_request = make_start_request(tmp_path)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start",
        MagicMock(side_effect=RuntimeError("backend failed")),
    )

    resp = core.deliver_sync_flush()

    assert resp.success is False
    assert resp.message == "Failed to prepare sync experiment"


def test_deliver_sync_flush_returns_failure_when_metrics_init_fails(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._start_record = make_start_record()
    core._start_request = make_start_request(tmp_path)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start",
        MagicMock(return_value=make_prepare_result(new_experiment=True)),
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.RunMetrics.new",
        MagicMock(side_effect=RuntimeError("metrics failed")),
    )

    resp = core.deliver_sync_flush()

    assert resp.success is False
    assert resp.message == "Failed to initialize sync metrics"


def test_read_uploads_newer_logs_and_captures_finish_record(tmp_path: Path):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    finish_record = FinishRecord(state=RunState.RUN_STATE_FINISHED, finished_at=make_timestamp())
    old_log = Record(log=LogRecord(epoch=1, level=LogLevel.LOG_LEVEL_INFO, line="old"))
    new_log = Record(log=LogRecord(epoch=3, level=LogLevel.LOG_LEVEL_INFO, line="new"))
    core._epoch = 1
    core._reader = FakeReader([Record(finish=finish_record), old_log, new_log])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert core._finish_record is not None
    assert core._finish_record.state == RunState.RUN_STATE_FINISHED
    assert [record.log.line for record in transport.records] == ["new"]
    assert core._epoch == 3


def test_read_skips_columns_already_in_metrics(tmp_path: Path):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    core._reader = FakeReader(  # type: ignore[assignment]
        [Record(column=ColumnRecord(column_key="loss", column_type=ColumnType.COLUMN_TYPE_SCALAR))]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    metrics = FakeMetrics()
    metrics.metrics["loss"] = FakeMetric()
    core._metrics = metrics  # type: ignore[assignment]

    asyncio.run(core.read())

    assert transport.records == []


def test_read_auto_defines_scalar_column_when_metric_is_missing(tmp_path: Path):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    scalar = ScalarRecord(
        key="loss",
        step=1,
        type=ColumnType.COLUMN_TYPE_SCALAR,
        value=ScalarValue(number=0.3),
    )
    core._reader = FakeReader([Record(scalar=scalar)])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    metrics = FakeMetrics()
    core._metrics = metrics  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record_kind(record) for record in transport.records] == ["column", "scalar"]
    assert transport.records[0].column.column_key == "loss"
    assert transport.records[0].column.column_type == ColumnType.COLUMN_TYPE_SCALAR
    assert metrics.get("loss") is not None


def test_read_auto_defines_media_column_when_metric_is_missing(tmp_path: Path):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    media = MediaRecord(key="images/sample", step=1, type=ColumnType.COLUMN_TYPE_IMAGE)
    core._reader = FakeReader([Record(media=media)])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    metrics = FakeMetrics()
    core._metrics = metrics  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record_kind(record) for record in transport.records] == ["column", "media"]
    assert transport.records[0].column.column_key == "images/sample"
    assert transport.records[0].column.column_type == ColumnType.COLUMN_TYPE_IMAGE
    assert metrics.get("images/sample") is not None


def test_read_skips_invalid_protobuf_payload_and_continues(tmp_path: Path):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    first = Record(
        scalar=ScalarRecord(
            key="loss",
            step=1,
            type=ColumnType.COLUMN_TYPE_SCALAR,
            value=ScalarValue(number=0.3),
        )
    )
    second = Record(
        scalar=ScalarRecord(
            key="loss",
            step=2,
            type=ColumnType.COLUMN_TYPE_SCALAR,
            value=ScalarValue(number=0.2),
        )
    )
    core._reader = FakeRawReader([first.SerializeToString(), b"\x80", second.SerializeToString()])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record_kind(record) for record in transport.records] == ["column", "scalar", "scalar"]
    assert [(record.scalar.key, record.scalar.step) for record in transport.records if record.HasField("scalar")] == [
        ("loss", 1),
        ("loss", 2),
    ]


def test_confirm_sync_finish_marks_missing_finish_record_as_crashed(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._read_executor = FakeExecutor()  # type: ignore[assignment]
    core._reader = FakeReader([])  # type: ignore[assignment]
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is True
    stop_mock.assert_called_once()
    call_kwargs = stop_mock.call_args
    assert call_kwargs.kwargs["state"] == RunState.RUN_STATE_CRASHED
    assert transport.finished is True
    assert [record_kind(record) for record in transport.records] == ["log"]
    assert "before writing a finish record" in transport.records[0].log.line


def test_confirm_sync_finish_uploads_error_log_for_failed_finish_record(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._read_executor = FakeExecutor()  # type: ignore[assignment]
    core._reader = FakeReader([])  # type: ignore[assignment]
    finished_at = make_timestamp()
    core._finish_record = FinishRecord(
        state=RunState.RUN_STATE_ABORTED,
        error="user aborted run",
        finished_at=finished_at,
    )
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is True
    stop_mock.assert_called_once()
    call_kwargs = stop_mock.call_args
    assert call_kwargs.kwargs["state"] == RunState.RUN_STATE_ABORTED
    assert transport.finished is True
    assert [record_kind(record) for record in transport.records] == ["log"]
    assert transport.records[0].log.level == LogLevel.LOG_LEVEL_ERROR
    assert transport.records[0].log.line == "user aborted run"


def test_confirm_sync_finish_uses_existing_finish_record(tmp_path: Path, monkeypatch):
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._read_executor = FakeExecutor()  # type: ignore[assignment]
    core._reader = FakeReader([])  # type: ignore[assignment]
    finished_at = make_timestamp()
    core._finish_record = FinishRecord(state=RunState.RUN_STATE_FINISHED, finished_at=finished_at)
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is True
    stop_mock.assert_called_once()
    call_kwargs = stop_mock.call_args
    assert call_kwargs.kwargs["state"] == RunState.RUN_STATE_FINISHED
    assert call_kwargs.kwargs["finished_at"] == finished_at
    assert transport.records == []
    assert transport.finished is True
