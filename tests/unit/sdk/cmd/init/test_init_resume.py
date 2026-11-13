"""
@file: test_init_resume.py
@description: 测试显式 run.dir 复用与跨进程续传前置校验语义。
"""

from pathlib import Path
from typing import List, Optional

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState, OperationStats
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import StartRecord
from swanlab.sdk.cmd import init as init_module
from swanlab.sdk.cmd.init import ensure_run_dir, init
from swanlab.sdk.internal.core_python.store import DataStoreWriter
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.protocol.core import NullCoreProtocol
from swanlab.sdk.typings.run import ModeType


def _write_run_file(
    run_dir: Path,
    run_id: str,
    *,
    filename: Optional[str] = None,
    with_start: bool = True,
    extra_records: int = 0,
) -> Path:
    run_file = run_dir / (filename or f"run-{run_id}.swanlab")
    writer = DataStoreWriter()
    writer.open(run_file)
    try:
        if with_start:
            rec = Record(num=-1, start=StartRecord(id=run_id, project="test-project"))
            writer.write(rec.SerializeToString())
        for i in range(extra_records):
            writer.write(Record(num=i + 1).SerializeToString())
    finally:
        writer.close()
    return run_file


class _FakeResumeCore(NullCoreProtocol):
    """可控模拟 Core：记录 DeliverRunStart 请求，模拟续传步数或已完成拒绝。"""

    def __init__(
        self,
        mode: ModeType,
        *,
        start_error: Optional[str] = None,
        global_step: int = -1,
        global_system_step: int = -1,
        new_experiment: bool = True,
        started_requests: Optional[List[DeliverRunStartRequest]] = None,
    ) -> None:
        super().__init__(mode)
        self._start_error = start_error
        self._global_step = global_step
        self._global_system_step = global_system_step
        self._new_experiment = new_experiment
        self._started_requests = started_requests if started_requests is not None else []

    def _start_any(self, req: DeliverRunStartRequest) -> DeliverRunStartResponse:
        self._started_requests.append(req)
        if self._start_error is not None:
            return DeliverRunStartResponse(success=False, message=self._start_error)
        return DeliverRunStartResponse(
            success=True,
            message="OK",
            path="/user/proj/exp" if self._mode == "online" else "",
            run=req.start_record,
            name=req.start_record.name or "exp",
            global_step=self._global_step,
            global_system_step=self._global_system_step,
            new_experiment=self._new_experiment,
            session_id="sess-1",
        )

    def _start_when_local(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._start_any(start_request)

    def _start_when_offline(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._start_any(start_request)

    def _start_when_online(self, start_request: DeliverRunStartRequest) -> DeliverRunStartResponse:
        return self._start_any(start_request)

    def _finish_when_local(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        return DeliverRunFinishResponse(success=True, message="OK")

    def _finish_when_offline(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        return DeliverRunFinishResponse(success=True, message="OK")

    def _finish_when_online(self, finish_request: DeliverRunFinishRequest) -> DeliverRunFinishResponse:
        return DeliverRunFinishResponse(success=True, message="OK")

    def _get_operation_when_online(self) -> GetOperationStatsResponse:
        return GetOperationStatsResponse(success=True, stats=OperationStats(state=CoreState.CORE_STATE_FINISHED))

    def _confirm_finish_when_enabled(self) -> ConfirmRunFinishResponse:
        return ConfirmRunFinishResponse(success=True, message="OK")


def _install_fake_core(
    monkeypatch,
    *,
    start_error: Optional[str] = None,
    global_step: int = -1,
    global_system_step: int = -1,
    new_experiment: bool = True,
) -> List[DeliverRunStartRequest]:
    requests_seen: List[DeliverRunStartRequest] = []

    def fake_create_core(ctx) -> _FakeResumeCore:
        return _FakeResumeCore(
            ctx.config.settings.mode,
            start_error=start_error,
            global_step=global_step,
            global_system_step=global_system_step,
            new_experiment=new_experiment,
            started_requests=requests_seen,
        )

    monkeypatch.setattr(init_module, "create_core", fake_create_core)
    monkeypatch.setattr("swanlab.sdk.internal.context.create_core", fake_create_core)
    return requests_seen


def _settings_with_dir(log_dir: Path, run_dir_name: Optional[str] = None) -> Settings:
    data: dict = {"log_dir": log_dir, "probe": {"monitor": False}}
    if run_dir_name is not None:
        data["run"] = {"dir": run_dir_name}
    return Settings.model_validate(data)


class TestEnsureRunDirReuse:
    def test_default_reuse_existing_false_still_raises(self, tmp_path):
        (tmp_path / "custom-dir").mkdir()
        with pytest.raises(FileExistsError):
            ensure_run_dir(tmp_path, "run-1", dir_name="custom-dir")

    def test_reuse_existing_true_keeps_existing_contents(self, tmp_path):
        target = tmp_path / "custom-dir"
        (target / "media").mkdir(parents=True)
        marker = target / "media" / "keep.txt"
        marker.write_text("keep")

        resolved = ensure_run_dir(tmp_path, "run-1", dir_name="custom-dir", reuse_existing=True)

        assert resolved == target
        assert marker.read_text() == "keep"


class TestExplicitRunDirReuse:
    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_creates_when_explicit_run_dir_missing(self, tmp_path, monkeypatch, mode):
        log_dir = tmp_path / "swanlog"
        requests_seen = _install_fake_core(monkeypatch)

        run = init(mode=mode, id="run-new", settings=_settings_with_dir(log_dir, "fixed-dir"))
        try:
            assert run.dir == log_dir / "fixed-dir"
            assert run.dir.is_dir()
            assert len(requests_seen) == 1
        finally:
            run.finish()

    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_reuses_empty_existing_dir_without_cleaning(self, tmp_path, monkeypatch, mode):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "fixed-dir"
        (existing / "media").mkdir(parents=True)
        (existing / "files").mkdir()
        (existing / "debug").mkdir()
        keep_file = existing / "media" / "artifact.bin"
        keep_file.write_bytes(b"payload")

        infos: List[str] = []
        monkeypatch.setattr(console, "info", lambda *args, **kwargs: infos.append(" ".join(str(a) for a in args)))
        requests_seen = _install_fake_core(monkeypatch)

        run = init(mode=mode, id="run-1", settings=_settings_with_dir(log_dir, "fixed-dir"))
        try:
            assert run.dir == existing
            assert keep_file.read_bytes() == b"payload"
            assert len(requests_seen) == 1
            assert any(str(existing) in msg for msg in infos)
        finally:
            run.finish()

    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_reuses_dir_with_matching_start_record_and_continues_step(self, tmp_path, monkeypatch, mode):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "resume-dir"
        (existing / "media").mkdir(parents=True)
        keep_file = existing / "media" / "img.png"
        keep_file.write_bytes(b"png")
        run_file = _write_run_file(existing, "exp-001", extra_records=2)

        requests_seen = _install_fake_core(monkeypatch, global_step=99, global_system_step=12, new_experiment=False)

        run = init(mode=mode, id="exp-001", settings=_settings_with_dir(log_dir, "resume-dir"))
        try:
            assert run.dir == existing
            assert run_file.exists()
            assert keep_file.read_bytes() == b"png"
            assert run._ctx.next_step() == 100
            assert len(requests_seen) == 1
            assert requests_seen[0].core_settings.run_dir == str(existing)
        finally:
            run.finish()

    def test_reuses_dir_with_env_var(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "env-dir"
        existing.mkdir(parents=True)
        _write_run_file(existing, "exp-env")

        monkeypatch.setenv("SWANLAB_RUN_DIR", "env-dir")
        requests_seen = _install_fake_core(monkeypatch, global_step=5, new_experiment=False)

        run = init(
            mode="offline", id="exp-env", log_dir=str(log_dir), settings=Settings(probe=Settings.Probe(monitor=False))
        )
        try:
            assert run.dir == existing
            assert run._ctx.next_step() == 6
            assert len(requests_seen) == 1
        finally:
            run.finish()

    def test_rejects_mismatched_run_id_before_starting_core(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "other-dir"
        existing.mkdir(parents=True)
        _write_run_file(existing, "exp-old")

        requests_seen = _install_fake_core(monkeypatch)

        with pytest.raises(RuntimeError, match="belongs to another run"):
            init(mode="offline", id="exp-new", settings=_settings_with_dir(log_dir, "other-dir"))

        assert requests_seen == []

    def test_propagates_core_already_finished_error(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "finished-dir"
        existing.mkdir(parents=True)
        _write_run_file(existing, "exp-done")

        _install_fake_core(
            monkeypatch,
            start_error="invalid local run file: run 'exp-done' is already finished, use `swanlab sync` instead",
        )

        with pytest.raises(RuntimeError, match="swanlab sync"):
            init(mode="offline", id="exp-done", settings=_settings_with_dir(log_dir, "finished-dir"))

    def test_rejects_multiple_run_files_before_starting_core(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "multi-dir"
        existing.mkdir(parents=True)
        _write_run_file(existing, "exp-1", filename="run-a.swanlab")
        _write_run_file(existing, "exp-1", filename="run-b.swanlab")

        requests_seen = _install_fake_core(monkeypatch)

        with pytest.raises(RuntimeError, match="Multiple run-\\*\\.swanlab files"):
            init(mode="offline", id="exp-1", settings=_settings_with_dir(log_dir, "multi-dir"))

        assert requests_seen == []

    def test_rejects_corrupt_run_file_before_starting_core(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "corrupt-dir"
        existing.mkdir(parents=True)
        (existing / "run-exp-bad.swanlab").write_bytes(b"")

        requests_seen = _install_fake_core(monkeypatch)

        with pytest.raises(RuntimeError, match="Failed to read existing run file"):
            init(mode="offline", id="exp-bad", settings=_settings_with_dir(log_dir, "corrupt-dir"))

        assert requests_seen == []

    def test_allows_header_only_run_file(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        existing = log_dir / "header-only-dir"
        existing.mkdir(parents=True)
        _write_run_file(existing, "exp-hdr", with_start=False)

        requests_seen = _install_fake_core(monkeypatch)

        run = init(mode="offline", id="exp-hdr", settings=_settings_with_dir(log_dir, "header-only-dir"))
        try:
            assert run.dir == existing
            assert len(requests_seen) == 1
        finally:
            run.finish()


def _mock_online_client(monkeypatch) -> None:
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.internal.run.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.internal.run.client.reset", lambda: None)


class TestResumeControlPlaneSeparation:
    def test_without_explicit_dir_never_reuses_even_with_resume_allow_and_same_id(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        _install_fake_core(monkeypatch)
        _mock_online_client(monkeypatch)

        run1 = init(mode="online", id="exp-same", resume="allow", settings=_settings_with_dir(log_dir))
        dir1 = run1.dir
        run1.finish()

        run2 = init(mode="online", id="exp-same", resume="allow", settings=_settings_with_dir(log_dir))
        dir2 = run2.dir
        run2.finish()

        assert dir1 != dir2

    def test_resume_true_with_custom_id_without_run_dir_creates_new_dir(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        _install_fake_core(monkeypatch)
        _mock_online_client(monkeypatch)

        run1 = init(mode="online", id="exp-named", resume=True, settings=_settings_with_dir(log_dir))
        dir1 = run1.dir
        run1.finish()

        run2 = init(mode="online", id="exp-named", resume=True, settings=_settings_with_dir(log_dir))
        dir2 = run2.dir
        run2.finish()

        assert dir1 != dir2

    def test_resume_true_with_custom_id_and_explicit_run_dir_reuses_dir(self, tmp_path, monkeypatch):
        log_dir = tmp_path / "swanlog"
        _install_fake_core(monkeypatch)
        _mock_online_client(monkeypatch)

        run1 = init(mode="online", id="exp-named", resume=True, settings=_settings_with_dir(log_dir, "shared-dir"))
        dir1 = run1.dir
        _write_run_file(dir1, "exp-named")
        run1.finish()

        run2 = init(mode="online", id="exp-named", resume=True, settings=_settings_with_dir(log_dir, "shared-dir"))
        dir2 = run2.dir
        run2.finish()

        assert dir1 == dir2 == log_dir / "shared-dir"

    def test_online_resume_never_with_id_still_rejected(self, tmp_path, monkeypatch):
        _mock_online_client(monkeypatch)
        with pytest.raises(AssertionError, match="Run id should not be provided when resume=never"):
            init(
                mode="online", id="exp-bad", resume="never", settings=_settings_with_dir(tmp_path / "swanlog", "dir-x")
            )

    @pytest.mark.parametrize(
        ("mode", "skip_store"),
        [("disabled", False), ("online", True)],
    )
    def test_disabled_and_skip_store_do_not_create_or_reuse_dir(self, tmp_path, monkeypatch, mode, skip_store):
        log_dir = tmp_path / "swanlog"
        target = log_dir / "no-mkdir-dir"
        _install_fake_core(monkeypatch)
        if mode == "online":
            _mock_online_client(monkeypatch)

        settings = Settings.model_validate(
            {
                "log_dir": log_dir,
                "run": {"dir": "no-mkdir-dir"},
                "core": {"skip_store": skip_store},
                "probe": {"monitor": False},
            }
        )
        run = init(mode=mode, settings=settings)
        try:
            assert run.dir == target
            assert not target.exists()
        finally:
            run.finish()
