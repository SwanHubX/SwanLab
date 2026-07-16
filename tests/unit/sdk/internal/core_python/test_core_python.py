"""
@author: cunyue
@file: test_core_python.py
@time: 2026/4/19
@description: CorePython 单元测试

测试分组：
  - TestCorePythonStart  : 各模式 deliver_run_start 行为
  - TestCorePythonPublish: 各模式 publish 行为
  - TestCorePythonFinish : 各模式 deliver_run_finish 行为
  - TestCorePythonGuard  : 防御逻辑
"""

from unittest.mock import MagicMock, patch

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings as CoreSettingsPb
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


def make_core_ctx(tmp_path) -> CoreContext:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    return CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=run_dir,
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
        )
    )


def make_start_record() -> StartRecord:
    return StartRecord(id="test-run-id", project="test-project")


def make_start_request(tmp_path, record: StartRecord) -> DeliverRunStartRequest:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    return DeliverRunStartRequest(
        core_settings=CoreSettingsPb(
            run_id="test-run-id",
            run_dir=str(run_dir),
            section_rule=0,
            record_batch=10000,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
        ),
        start_record=record,
    )


def set_online_params(core: CorePython) -> None:
    core._ctx.set_online_params(
        username="test-user",
        project="test-project",
        project_id="test-project-id",
        project_version=1,
        experiment_id="test-experiment-id",
    )


# ============================================================
# TestCorePythonStart
# ============================================================


class TestCorePythonStart:
    @pytest.mark.parametrize("mode", ["disabled", "local", "offline"])
    def test_non_online_modes(self, tmp_path, mode):
        core = CorePython(mode)
        core._ctx = make_core_ctx(tmp_path)
        req = make_start_request(tmp_path, make_start_record())
        resp = core.deliver_run_start(req)

        assert resp.success is True
        assert core._transport is None
        if mode == "disabled":
            assert core._store is None
        else:
            assert core._store is not None

    def test_online_mode(self, tmp_path, monkeypatch):
        core = CorePython("online")
        record = make_start_record()
        mock_deliver = MagicMock(return_value=DeliverRunStartResponse(success=True, message="OK", run=record))

        def _report_run_start_and_set_attrs(rec):
            """mock _report_run_start，同时设置 online 模式所需的属性。"""
            set_online_params(core)
            core._metrics = MagicMock()
            return mock_deliver(rec)

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)

        with patch("swanlab.sdk.internal.core_python.core.Heartbeat") as MockHeartbeat:
            resp = core.deliver_run_start(make_start_request(tmp_path, record))

        assert resp.success is True
        assert core._store is not None
        assert core._transport is not None
        mock_deliver.assert_called_once_with(record)
        MockHeartbeat.assert_called_once_with("test-experiment-id")
        MockHeartbeat.return_value.start.assert_called_once()

    def test_start_failure_includes_original_error(self, tmp_path, monkeypatch):
        core = CorePython("local")
        error_message = "Cannot create SwanLab data store file /root/run.swanlab: permission denied."

        def deny_store_open(self, filename):
            raise PermissionError(error_message)

        monkeypatch.setattr("swanlab.sdk.internal.core_python.store.DataStoreWriter.open", deny_store_open)

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert "Failed to start run" in resp.message
        assert error_message in resp.message
        assert core._active is False


# ============================================================
# TestCorePythonFinish
# ============================================================


class TestCorePythonFinish:
    @pytest.mark.parametrize("mode", ["disabled", "local", "offline"])
    def test_non_online_modes(self, tmp_path, mode):
        core = CorePython(mode)
        core._ctx = make_core_ctx(tmp_path)
        core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        resp = core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))

        assert resp.success is True
        assert core._store is None

    def test_online_closes_store_and_transport(self, tmp_path, monkeypatch):
        core = CorePython("online")

        mock_start = MagicMock(
            return_value=DeliverRunStartResponse(success=True, message="OK", run=make_start_record())
        )

        def _report_run_start_and_set_attrs(rec):
            set_online_params(core)
            core._metrics = MagicMock()
            return mock_start(rec)

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)

        with patch("swanlab.sdk.internal.core_python.core.Heartbeat") as MockHeartbeat:
            core.deliver_run_start(make_start_request(tmp_path, make_start_record()))
        mock_heartbeat = MockHeartbeat.return_value

        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", lambda *a, **kw: None)
        core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))

        assert core._store is None
        assert core._transport is not None
        assert core._transport.join(timeout=2)
        confirm_resp = core.confirm_run_finish()
        assert core._transport is None
        mock_heartbeat.stop.assert_called_once()
        assert isinstance(confirm_resp, ConfirmRunFinishResponse)
        assert confirm_resp.success is True

    def test_confirm_run_finish_waits_without_timeout_before_marking_finished(self, tmp_path, monkeypatch):
        core = CorePython("online")
        core._ctx = make_core_ctx(tmp_path)
        set_online_params(core)
        tracker = UploadTracker()
        tracker.set_state(CoreState.CORE_STATE_RUNNING)
        transport = MagicMock()
        transport.finish.return_value = True
        core._tracker = tracker
        core._transport = transport
        core._pending_online_finish_record = FinishRecord()
        stop_experiment_mock = MagicMock()
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", stop_experiment_mock)

        resp = core.confirm_run_finish()

        assert isinstance(resp, ConfirmRunFinishResponse)
        assert resp.success is True
        assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED
        transport.finish.assert_called_once_with(timeout=None)
        stop_experiment_mock.assert_called_once()

    def test_confirm_run_finish_preserves_running_transport(self, tmp_path, monkeypatch):
        core = CorePython("online")
        core._ctx = make_core_ctx(tmp_path)
        set_online_params(core)
        transport = MagicMock()
        transport.finish.return_value = False
        heartbeat = MagicMock()
        stop_experiment_mock = MagicMock()
        core._transport = transport
        core._heartbeat = heartbeat
        core._pending_online_finish_record = FinishRecord()
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", stop_experiment_mock)

        resp = core.confirm_run_finish()

        assert isinstance(resp, ConfirmRunFinishResponse)
        assert resp.success is False
        assert core._transport is transport
        transport.finish.assert_called_once_with(timeout=None)
        heartbeat.stop.assert_not_called()
        stop_experiment_mock.assert_not_called()


# ============================================================
# TestCorePythonGuard
# ============================================================


class TestCorePythonGuard:
    def test_double_start_raises(self, tmp_path):
        core = CorePython("local")
        core._ctx = make_core_ctx(tmp_path)
        core.deliver_run_start(make_start_request(tmp_path, make_start_record()))
        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))
        assert resp.success is False
        assert "Run has already been active" in resp.message

    def test_finish_before_start(self, tmp_path):
        core = CorePython("local")
        core._ctx = make_core_ctx(tmp_path)

        resp = core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))
        assert resp.success is False
        assert "Run is not active" in resp.message

    def test_fork_raises(self, tmp_path):
        core = CorePython("disabled")

        with pytest.raises(RuntimeError, match="should not be called"):
            core.fork()
