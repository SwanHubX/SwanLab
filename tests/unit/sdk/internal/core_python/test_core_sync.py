import asyncio
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest, SyncRange
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord, ScalarValue
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, ResumeMode, RunState, StartRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel, LogRecord
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.metrics import RunMetrics
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.core_python.sync import CoreSyncPython, SyncRangeFilter
from swanlab.sdk.internal.core_python.utils import PrepareExperimentStartResult
from swanlab.utils.time import parse_timestamp_ns


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
        experiment={
            "cuid": "experiment-id",
            "slug": "sync-run-id",
            "name": "synced-run",
            "createdAt": "2024-08-01T00:00:00Z",
        },
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
        # 模拟上传线程异常退出：join 成功但缓冲区未排空
        self.drain_ok = True
        self.finish_error = "Upload did not finish: 1 record(s) retained (worker crashed: RuntimeError)."

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
        return self.join(timeout=timeout) and self.drain_ok

    def finish_summary(self) -> str:
        return self.finish_error


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
    def __init__(self, column_type=None, column=None):
        self.updated = []
        self.column_type = column_type
        self._column = column

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

    def allow_replay(self) -> None:
        pass

    def define_scalar(self, **kwargs) -> FakeMetric:
        column = kwargs.get("column")
        metric = FakeMetric(column_type=column.column_type if column else None, column=column)
        self.metrics[kwargs["key"]] = metric
        return metric

    def define_media(self, **kwargs) -> FakeMetric:
        column = kwargs.get("column")
        metric = FakeMetric(column_type=column.column_type if column else None, column=column)
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
    assert core._new_experiment is True
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
    get_summary.assert_called_once_with("project-id", "experiment-id", created_at="2024-08-01T00:00:00Z")


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
    assert resp.message.startswith("Failed to prepare sync experiment")
    assert "backend failed" in resp.message


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


def test_confirm_sync_finish_marks_new_experiment_as_crashed_without_finish_record(tmp_path: Path, monkeypatch):
    """新建实验且没有 FinishRecord：按 CRASHED 收尾，但不注入伪造的错误日志。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    core._new_experiment = True
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
    assert transport.records == []


def test_confirm_sync_finish_keeps_existing_experiment_running_without_finish_record(tmp_path: Path, monkeypatch):
    """已存在的实验且没有 FinishRecord：只结束本次同步，保持 RUNNING。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    core._new_experiment = False
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._read_executor = FakeExecutor()  # type: ignore[assignment]
    core._reader = FakeReader([])  # type: ignore[assignment]
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is True
    stop_mock.assert_not_called()
    assert transport.finished is True
    assert transport.records == []


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


# ============================================================
# define_metric: offline sync parity tests
# ============================================================


def test_sync_replays_column_with_xaxis_hidden(tmp_path: Path):
    """带 x_axis/hidden 的 ColumnRecord 在 sync 回放时保留字段。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    col = ColumnRecord(
        column_key="train/loss",
        column_type=ColumnType.COLUMN_TYPE_SCALAR,
        x_axis="train/epoch",
        section_name="Training",
        hidden=True,
    )
    core._reader = FakeReader([Record(column=col)])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert len(transport.records) == 1
    replayed = transport.records[0].column
    assert replayed.x_axis == "train/epoch"
    assert replayed.section_name == "Training"
    assert replayed.hidden is True


def test_sync_skips_when_key_already_exists(tmp_path: Path):
    """两条同 key column → 只上传第一条，第二条跳过（first-writer-wins）。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    col1 = ColumnRecord(column_key="loss", column_type=ColumnType.COLUMN_TYPE_SCALAR, section_name="A")
    col2 = ColumnRecord(column_key="loss", column_type=ColumnType.COLUMN_TYPE_SCALAR, section_name="B")
    core._reader = FakeReader([Record(column=col1), Record(column=col2)])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    col_records = [r.column for r in transport.records if record_kind(r) == "column"]
    assert len(col_records) == 1
    assert col_records[0].section_name == "A"


# ============================================================
# range sync: step / time filtering and in-place replay
# ============================================================


def make_timestamp_at(seconds: int, nanos: int = 0) -> Timestamp:
    ts = Timestamp()
    ts.seconds = seconds
    ts.nanos = nanos
    return ts


def make_range_core(tmp_path: Path, sync_range: SyncRange) -> CoreSyncPython:
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    core._filter = SyncRangeFilter.from_proto(sync_range)
    return core


def test_read_filters_scalars_by_step_range(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_step=2, end_step=2))
    records = [
        Record(scalar=ScalarRecord(key="loss", step=step, type=ColumnType.COLUMN_TYPE_SCALAR)) for step in (1, 2, 3)
    ]
    core._reader = FakeReader(records)  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.scalar.step for record in transport.records if record.HasField("scalar")] == [2]


def test_read_filters_scalars_by_time_range(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_time=make_timestamp_at(100), end_time=make_timestamp_at(200)))
    records = [
        Record(scalar=ScalarRecord(key="loss", step=step, type=ColumnType.COLUMN_TYPE_SCALAR, timestamp=ts))
        for step, ts in (
            (1, make_timestamp_at(99)),
            (2, make_timestamp_at(100)),
            (3, make_timestamp_at(200)),
            (4, make_timestamp_at(201)),
        )
    ]
    core._reader = FakeReader(records)  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.scalar.step for record in transport.records if record.HasField("scalar")] == [2, 3]


def test_read_falls_back_to_record_timestamp(tmp_path: Path):
    """指标自身没有时间戳时，使用记录落盘时间参与时间过滤。"""
    core = make_range_core(tmp_path, SyncRange(start_time=make_timestamp_at(100)))
    core._reader = FakeReader(  # type: ignore[assignment]
        [
            Record(
                timestamp=make_timestamp_at(150),
                scalar=ScalarRecord(key="loss", step=1, type=ColumnType.COLUMN_TYPE_SCALAR),
            ),
            Record(
                timestamp=make_timestamp_at(50),
                scalar=ScalarRecord(key="loss", step=2, type=ColumnType.COLUMN_TYPE_SCALAR),
            ),
        ]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.scalar.step for record in transport.records if record.HasField("scalar")] == [1]


def test_read_filters_media_by_step_range(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(end_step=1))
    media = MediaRecord(key="images/sample", step=2, type=ColumnType.COLUMN_TYPE_IMAGE)
    core._reader = FakeReader([Record(media=media)])  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert transport.records == []


def test_read_filters_logs_by_time_range(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_time=make_timestamp_at(100), end_time=make_timestamp_at(200)))
    records = [
        Record(log=LogRecord(epoch=1, level=LogLevel.LOG_LEVEL_INFO, line="old", timestamp=make_timestamp_at(99))),
        Record(log=LogRecord(epoch=1, level=LogLevel.LOG_LEVEL_INFO, line="in", timestamp=make_timestamp_at(150))),
        Record(log=LogRecord(epoch=2, level=LogLevel.LOG_LEVEL_INFO, line="new", timestamp=make_timestamp_at(201))),
    ]
    core._reader = FakeReader(records)  # type: ignore[assignment]
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.log.line for record in transport.records] == ["in"]
    assert core._epoch == 1


def test_read_skips_logs_when_only_step_range_is_set(tmp_path: Path):
    """step 范围只约束指标，未设置时间范围时不回放终端日志。"""
    core = make_range_core(tmp_path, SyncRange(start_step=1))
    core._reader = FakeReader(  # type: ignore[assignment]
        [Record(log=LogRecord(epoch=1, level=LogLevel.LOG_LEVEL_INFO, line="line"))]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert transport.records == []


def test_read_keeps_old_log_epochs_in_range_mode(tmp_path: Path):
    """范围回放不受云端已有 epoch 下界限制，保留原始 epoch 以便原地恢复对齐。"""
    core = make_range_core(tmp_path, SyncRange(start_time=make_timestamp_at(100)))
    core._epoch = 10
    core._reader = FakeReader(  # type: ignore[assignment]
        [
            Record(
                log=LogRecord(epoch=3, level=LogLevel.LOG_LEVEL_INFO, line="old line", timestamp=make_timestamp_at(150))
            )
        ]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.log.epoch for record in transport.records] == [3]
    assert core._epoch == 10


def test_read_skips_save_records_in_inplace_replay(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_step=0))
    core._in_place = True
    core._reader = FakeReader(  # type: ignore[assignment]
        [Record(save=SaveRecord(name="config", type=SaveType.SAVE_TYPE_CONFIG))]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert transport.records == []


def test_read_keeps_save_records_when_range_creates_new_experiment(tmp_path: Path):
    """范围同步新建实验时仍需要 config/metadata，否则实验缺少元信息。"""
    core = make_range_core(tmp_path, SyncRange(start_step=0))
    core._in_place = False
    core._reader = FakeReader(  # type: ignore[assignment]
        [Record(save=SaveRecord(name="config", type=SaveType.SAVE_TYPE_CONFIG))]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record_kind(record) for record in transport.records] == ["save"]


def test_read_keeps_columns_in_range_mode(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_step=100))
    core._reader = FakeReader(  # type: ignore[assignment]
        [Record(column=ColumnRecord(column_key="loss", column_type=ColumnType.COLUMN_TYPE_SCALAR))]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._metrics = FakeMetrics()  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record_kind(record) for record in transport.records] == ["column"]


def test_read_replays_steps_already_in_remote_summary(tmp_path: Path):
    """范围同步重置摘要中的 min_step，允许回放已存在的 step。"""
    core = make_range_core(tmp_path, SyncRange(start_step=1, end_step=3))
    metrics, _, _, _ = RunMetrics.new(
        {"log": None, "media": None, "scalar": [{"key": "loss", "step": 5}]}, ctx=core._ctx
    )
    metrics.allow_replay()
    core._metrics = metrics
    core._reader = FakeReader(  # type: ignore[assignment]
        [
            Record(
                scalar=ScalarRecord(
                    key="loss", step=3, type=ColumnType.COLUMN_TYPE_SCALAR, value=ScalarValue(number=1.0)
                )
            )
        ]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record.scalar.step for record in transport.records if record.HasField("scalar")] == [3]


def test_read_skips_step_at_or_below_summary_min_step_without_replay(tmp_path: Path):
    """未开启范围同步时，摘要中的 min_step 仍然拦截旧 step。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    metrics, _, _, _ = RunMetrics.new(
        {"log": None, "media": None, "scalar": [{"key": "loss", "step": 5}]}, ctx=core._ctx
    )
    core._metrics = metrics
    core._reader = FakeReader(  # type: ignore[assignment]
        [
            Record(
                scalar=ScalarRecord(
                    key="loss", step=3, type=ColumnType.COLUMN_TYPE_SCALAR, value=ScalarValue(number=1.0)
                )
            )
        ]
    )
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]

    asyncio.run(core.read())

    assert [record for record in transport.records if record.HasField("scalar")] == []


def test_deliver_sync_start_rejects_invalid_step_range(tmp_path: Path):
    write_run_file(tmp_path, Record(start=make_start_record()))
    core = CoreSyncPython()

    resp = core.deliver_sync_start(make_start_request(tmp_path, range=SyncRange(start_step=10, end_step=5)))

    assert resp.success is False
    assert "start_step is greater than end_step" in resp.message


def test_deliver_sync_start_rejects_invalid_time_range(tmp_path: Path):
    write_run_file(tmp_path, Record(start=make_start_record()))
    core = CoreSyncPython()

    resp = core.deliver_sync_start(
        make_start_request(
            tmp_path,
            range=SyncRange(start_time=make_timestamp_at(200), end_time=make_timestamp_at(100)),
        )
    )

    assert resp.success is False
    assert "start_time is later than end_time" in resp.message


def test_queue_finish_records_skipped_in_range_mode(tmp_path: Path):
    core = make_range_core(tmp_path, SyncRange(start_step=1))
    transport = FakeTransport(core._ctx)
    core._transport = transport  # type: ignore[assignment]
    core._finish_record = FinishRecord(
        state=RunState.RUN_STATE_CRASHED, error="boom", finished_at=make_timestamp_at(100)
    )

    core._queue_finish_records()

    assert transport.records == []


def make_confirm_core(tmp_path: Path) -> CoreSyncPython:
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._ctx.set_online_params("alice", "demo", "project-id", None, "experiment-id")
    core._read_executor = FakeExecutor()  # type: ignore[assignment]
    core._reader = FakeReader([])  # type: ignore[assignment]
    return core


def test_confirm_sync_finish_reports_upload_failure(tmp_path: Path, monkeypatch):
    """上传线程异常退出时，全量同步必须报告失败，而不是当作成功结束实验。"""
    core = make_confirm_core(tmp_path)
    transport = FakeTransport(core._ctx)
    transport.drain_ok = False
    core._transport = transport  # type: ignore[assignment]
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is False
    assert "retained" in resp.message
    stop_mock.assert_not_called()


def test_confirm_sync_finish_in_place_reports_upload_failure(tmp_path: Path, monkeypatch):
    """原地回填同样要报告排空失败，不能用成功掩盖未上传的数据。"""
    core = make_confirm_core(tmp_path)
    core._in_place = True
    transport = FakeTransport(core._ctx)
    transport.drain_ok = False
    core._transport = transport  # type: ignore[assignment]
    stop_mock = MagicMock()
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", stop_mock)

    resp = core.confirm_sync_finish()

    assert resp.success is False
    assert "retained" in resp.message
    stop_mock.assert_not_called()


def test_deliver_sync_flush_range_resumes_existing_experiment(tmp_path: Path, monkeypatch):
    """范围同步也统一走创建/恢复流程，已有实验由 Server 执行 resume。"""
    core = CoreSyncPython()
    core._ctx = CoreContext(config=make_core_config(tmp_path), mode="sync")
    core._start_record = make_start_record()
    core._start_request = make_start_request(tmp_path, range=SyncRange(start_step=1))
    core._filter = SyncRangeFilter.from_proto(SyncRange(start_step=1))
    prepare = MagicMock(return_value=make_prepare_result(new_experiment=False))
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.prepare_experiment_start", prepare)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.get_experiment_summary",
        MagicMock(return_value={"log": None, "media": None, "scalar": None}),
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.RunMetrics.new", MagicMock(return_value=(FakeMetrics(), -1, -1, -1))
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.Transport",
        lambda ctx, tracker=None: FakeTransport(ctx, tracker),
    )
    core._read_executor = FakeExecutor()  # type: ignore[assignment]

    assert core.deliver_sync_flush().success is True
    prepare.assert_called_once()
    assert core._new_experiment is False
    assert core._in_place is True


# ============================================================
# range filter: step / time boundary matching
# ============================================================


def test_filter_disabled_without_range():
    f = SyncRangeFilter.from_proto(None)

    assert f.enabled is False
    assert f.match(step=1, timestamp=make_timestamp_at(1)) is True


def test_filter_disabled_when_timestamps_are_zero():
    f = SyncRangeFilter.from_proto(SyncRange(start_time=Timestamp(), end_time=Timestamp()))

    assert f.enabled is False


def test_filter_step_range_is_inclusive():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=10, end_step=20))

    assert f.enabled is True
    assert f.match(step=10) is True
    assert f.match(step=20) is True
    assert f.match(step=9) is False
    assert f.match(step=21) is False


def test_filter_step_range_ignores_records_without_step():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=10, end_step=20))

    # 终端日志没有 step，只受时间范围约束
    assert f.match(timestamp=make_timestamp_at(1)) is True


def test_filter_reports_time_bounds():
    assert SyncRangeFilter.from_proto(SyncRange(start_step=1)).has_time is False
    assert SyncRangeFilter.from_proto(SyncRange(start_time=make_timestamp_at(1))).has_time is True


def test_filter_open_ended_step_range():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=100))

    assert f.match(step=100) is True
    assert f.match(step=10_000) is True
    assert f.match(step=99) is False


def test_filter_time_range_is_inclusive():
    f = SyncRangeFilter.from_proto(SyncRange(start_time=make_timestamp_at(100, 500), end_time=make_timestamp_at(200)))

    assert f.match(step=1, timestamp=make_timestamp_at(100, 500)) is True
    assert f.match(step=1, timestamp=make_timestamp_at(200)) is True
    assert f.match(step=1, timestamp=make_timestamp_at(100, 499)) is False
    assert f.match(step=1, timestamp=make_timestamp_at(201)) is False


def test_filter_combines_step_and_time():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=2, end_step=3, start_time=make_timestamp_at(100)))

    assert f.match(step=2, timestamp=make_timestamp_at(100)) is True
    assert f.match(step=1, timestamp=make_timestamp_at(100)) is False
    assert f.match(step=2, timestamp=make_timestamp_at(99)) is False


def test_filter_keeps_subsecond_endpoint():
    """从同一字面量解析出的时间戳必须命中它自己构成的闭区间。"""
    instant = parse_timestamp_ns("2024-08-01T00:00:00.000001Z")
    timestamp = Timestamp()
    timestamp.FromNanoseconds(instant)

    f = SyncRangeFilter(start_ns=instant, end_ns=instant)

    assert f.match(timestamp=timestamp) is True
