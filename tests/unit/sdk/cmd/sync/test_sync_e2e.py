from pathlib import Path
from typing import Iterable

import pytest

import swanlab
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord, ScalarValue
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import RunState, StartRecord
from swanlab.sdk.cmd import sync as sync_cmd
from swanlab.sdk.internal.core_python.store import LEVELDBLOG_HEADER_LEN, DataStoreReader, DataStoreWriter
from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.core_python.utils import PrepareExperimentStartResult
from swanlab.sdk.internal.settings import Settings


def read_records(run_dir: Path) -> list[Record]:
    files = sorted(run_dir.glob("run-*.swanlab"))
    assert len(files) == 1

    reader = DataStoreReader()
    reader.open(files[0])
    records = []
    for record_bytes in reader:
        record = Record()
        record.ParseFromString(record_bytes)
        records.append(record)
    reader.close()
    return records


def record_kinds(records: Iterable[Record]) -> list[str]:
    kinds = []
    for record in records:
        for kind in [
            "start",
            "finish",
            "column",
            "scalar",
            "media",
            "config",
            "log",
            "metadata",
            "requirements",
            "conda",
            "save",
        ]:
            if record.HasField(kind):
                kinds.append(kind)
                break
    return kinds


def scalar_signatures(records: Iterable[Record]) -> set[tuple[str, int]]:
    return {(record.scalar.key, record.scalar.step) for record in records if record.HasField("scalar")}


def column_signatures(records: Iterable[Record]) -> set[str]:
    return {record.column.column_key for record in records if record.HasField("column")}


def make_prepare_result(new_experiment: bool = True) -> PrepareExperimentStartResult:
    return PrepareExperimentStartResult(
        username="alice",
        project="demo",
        project_info={
            "cuid": "project-id",
            "name": "demo",
            "username": "alice",
            "path": "/alice/demo",
            "visibility": "PRIVATE",
            "_count": {"experiments": 0, "contributors": 0, "collaborators": 0, "clones": 0},
        },
        experiment={"cuid": "experiment-id", "slug": "sync-e2e-run", "name": "sync-e2e-run"},
        new_experiment=new_experiment,
        name="sync-e2e-run",
        color="#abcdef",
    )


def create_offline_run() -> Path:
    run = swanlab.init(mode="offline", project="sync-e2e", id="sync-e2e-run")
    run_dir = run._ctx.run_dir
    swanlab.log({"loss": 0.2, "acc": 0.8}, step=1)
    swanlab.log({"loss": 0.1, "acc": 0.9}, step=2)
    swanlab.finish()
    return run_dir


def make_start_record() -> StartRecord:
    return StartRecord(id="sync-e2e-run", workspace="alice", project="demo")


def make_scalar_record(step: int, value: float) -> Record:
    return Record(
        scalar=ScalarRecord(
            key="loss",
            step=step,
            type=ColumnType.COLUMN_TYPE_SCALAR,
            value=ScalarValue(number=value),
        )
    )


def write_raw_run(run_dir: Path, records: list[Record]) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    run_file = run_dir / "run-sync-e2e-run.swanlab"
    writer = DataStoreWriter()
    writer.open(run_file)
    for record in records:
        writer.write(record.SerializeToString())
    writer.close()
    return run_file


def corrupt_record_payload(run_file: Path, payloads: list[bytes], record_index: int) -> None:
    offset = LEVELDBLOG_HEADER_LEN
    for payload in payloads[:record_index]:
        offset += LEVELDBLOG_HEADER_LEN + len(payload)
    offset += LEVELDBLOG_HEADER_LEN

    data = bytearray(run_file.read_bytes())
    data[offset] ^= 0xFF
    run_file.write_bytes(bytes(data))


class FakeTransport:
    def __init__(self, ctx):
        self.ctx = ctx
        self.records: list[Record] = []
        self.started = False
        self.finished = False

    def put(self, records: list[Record]) -> None:
        self.records.extend(records)

    def start(self) -> None:
        self.started = True

    def finish(self) -> None:
        self.finished = True


class ImmediateExecutor:
    def __init__(self):
        self.started = False
        self.waited = False
        self.closed = False

    def start(self, coro) -> None:
        self.started = True
        try:
            coro.send(None)
        except StopIteration:
            pass

    def wait(self) -> None:
        self.waited = True

    def close(self) -> None:
        self.closed = True


def test_e2e_offline_run_writes_expected_records_to_swanlab_file():
    run_dir = create_offline_run()

    records = read_records(run_dir)
    kinds = record_kinds(records)

    assert "start" in kinds
    assert "column" in kinds
    assert "scalar" in kinds
    assert "finish" in kinds
    assert {("loss", 1), ("loss", 2), ("acc", 1), ("acc", 2)} <= scalar_signatures(records)
    assert {"loss", "acc"} <= column_signatures(records)


def test_e2e_sync_uploads_records_from_generated_swanlab_file(monkeypatch):
    run_dir = create_offline_run()
    source_records = read_records(run_dir)
    transports: list[FakeTransport] = []
    stopped = []

    core = CoreSyncPython()
    core._read_executor = ImmediateExecutor()  # type: ignore[assignment]

    def make_transport(ctx) -> FakeTransport:
        transport = FakeTransport(ctx)
        transports.append(transport)
        return transport

    def fake_stop_experiment(username, project, experiment_id, *, state, finished_at):
        stopped.append((username, project, experiment_id, state, finished_at))

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start", lambda record: make_prepare_result()
    )
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.Transport", make_transport)
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", fake_stop_experiment)

    sync_cmd.sync(
        run_dir,
        settings=Settings(
            api_key="test-api-key",
            project=Settings.Project(workspace="alice", name="demo"),
            run=Settings.Run(id="sync-e2e-run"),
        ),
    )

    assert len(transports) == 1
    uploaded_records = transports[0].records
    assert transports[0].started is True
    assert transports[0].finished is True
    assert {("loss", 1), ("loss", 2), ("acc", 1), ("acc", 2)} <= scalar_signatures(uploaded_records)
    assert {"loss", "acc"} <= column_signatures(uploaded_records)
    assert "start" not in record_kinds(uploaded_records)
    assert "finish" not in record_kinds(uploaded_records)
    assert scalar_signatures(source_records) <= scalar_signatures(uploaded_records)
    assert stopped
    assert stopped[0][0:3] == ("alice", "demo", "experiment-id")


def test_e2e_sync_skips_records_already_in_remote_summary(monkeypatch):
    run_dir = create_offline_run()
    transports: list[FakeTransport] = []

    core = CoreSyncPython()
    core._read_executor = ImmediateExecutor()  # type: ignore[assignment]

    def make_transport(ctx) -> FakeTransport:
        transport = FakeTransport(ctx)
        transports.append(transport)
        return transport

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start",
        lambda record: make_prepare_result(new_experiment=False),
    )
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.get_experiment_summary",
        lambda project_id, experiment_id: {
            "log": None,
            "media": None,
            "scalar": [{"key": "loss", "step": 1}, {"key": "acc", "step": 1}],
        },
    )
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.Transport", make_transport)
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", lambda *args, **kwargs: None)

    sync_cmd.sync(
        run_dir,
        settings=Settings(
            api_key="test-api-key",
            project=Settings.Project(workspace="alice", name="demo"),
            run=Settings.Run(id="sync-e2e-run"),
        ),
    )

    assert len(transports) == 1
    assert scalar_signatures(transports[0].records) == {("loss", 2), ("acc", 2)}


def test_e2e_sync_stops_uploading_after_corrupted_record(tmp_path: Path, monkeypatch):
    records = [
        Record(start=make_start_record()),
        make_scalar_record(step=1, value=0.3),
        make_scalar_record(step=2, value=0.2),
        make_scalar_record(step=3, value=0.1),
    ]
    payloads = [record.SerializeToString() for record in records]
    run_file = write_raw_run(tmp_path, records)
    corrupt_record_payload(run_file, payloads, record_index=2)
    transports: list[FakeTransport] = []

    core = CoreSyncPython()
    core._read_executor = ImmediateExecutor()  # type: ignore[assignment]

    def make_transport(ctx) -> FakeTransport:
        transport = FakeTransport(ctx)
        transports.append(transport)
        return transport

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start", lambda record: make_prepare_result()
    )
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.Transport", make_transport)
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", lambda *args, **kwargs: None)

    sync_cmd.sync(
        tmp_path,
        settings=Settings(
            api_key="test-api-key",
            project=Settings.Project(workspace="alice", name="demo"),
            run=Settings.Run(id="sync-e2e-run"),
        ),
    )

    assert len(transports) == 1
    assert scalar_signatures(transports[0].records) == {("loss", 1)}


def test_e2e_sync_marks_missing_finish_record_as_crashed(tmp_path: Path, monkeypatch):
    write_raw_run(tmp_path, [Record(start=make_start_record()), make_scalar_record(step=1, value=0.3)])
    transports: list[FakeTransport] = []
    stopped = []

    core = CoreSyncPython()
    core._read_executor = ImmediateExecutor()  # type: ignore[assignment]

    def make_transport(ctx) -> FakeTransport:
        transport = FakeTransport(ctx)
        transports.append(transport)
        return transport

    def fake_stop_experiment(username, project, experiment_id, *, state, finished_at):
        stopped.append((username, project, experiment_id, state, finished_at))

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.internal.core_python.sync.prepare_experiment_start", lambda record: make_prepare_result()
    )
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.Transport", make_transport)
    monkeypatch.setattr("swanlab.sdk.internal.core_python.sync.stop_experiment", fake_stop_experiment)

    sync_cmd.sync(
        tmp_path,
        settings=Settings(
            api_key="test-api-key",
            project=Settings.Project(workspace="alice", name="demo"),
            run=Settings.Run(id="sync-e2e-run"),
        ),
    )

    assert len(transports) == 1
    assert stopped[0][0:4] == ("alice", "demo", "experiment-id", RunState.RUN_STATE_CRASHED)
    assert "log" in record_kinds(transports[0].records)
    error_logs = [record.log.line for record in transports[0].records if record.HasField("log")]
    assert any("before writing a finish record" in line for line in error_logs)


def test_e2e_sync_fails_with_legacy_version_file(tmp_path: Path, monkeypatch):
    """sync 遇到旧版本（v1）文件时应失败并给出版本错误提示。"""
    import struct

    from swanlab.sdk.internal.core_python.store import LEVELDBLOG_HEADER_IDENT, LEVELDBLOG_HEADER_MAGIC

    run_file = tmp_path / "run-legacy-run.swanlab"
    with open(run_file, "wb") as f:
        f.write(struct.pack("<4sHB", LEVELDBLOG_HEADER_IDENT, LEVELDBLOG_HEADER_MAGIC, 1))

    core = CoreSyncPython()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="version"):
        sync_cmd.sync(
            tmp_path,
            settings=Settings(
                api_key="test-api-key",
                project=Settings.Project(workspace="alice", name="demo"),
                run=Settings.Run(id="legacy-run"),
            ),
        )
