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
    DeliverRunFinishRequest,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
)
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings as CoreSettingsPb
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext


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
            core._username = "test-user"
            core._project = "test-project"
            core._project_id = "test-project-id"
            core._experiment_id = "test-experiment-id"
            core._metrics = MagicMock()
            return mock_deliver(rec)

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)

        with patch("swanlab.sdk.internal.core_python.Heartbeat") as MockHeartbeat:
            resp = core.deliver_run_start(make_start_request(tmp_path, record))

        assert resp.success is True
        assert core._store is not None
        assert core._transport is not None
        mock_deliver.assert_called_once_with(record)
        MockHeartbeat.assert_called_once_with("test-experiment-id")
        MockHeartbeat.return_value.start.assert_called_once()


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
            """mock _report_run_start，同时设置 online 模式所需的属性。"""
            core._username = "test-user"
            core._project = "test-project"
            core._project_id = "test-project-id"
            core._experiment_id = "test-experiment-id"
            core._metrics = MagicMock()
            return mock_start(rec)

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)

        with patch("swanlab.sdk.internal.core_python.Heartbeat") as MockHeartbeat:
            core.deliver_run_start(make_start_request(tmp_path, make_start_record()))
        mock_heartbeat = MockHeartbeat.return_value

        mock_finish = MagicMock(return_value=None)
        monkeypatch.setattr(core, "_report_run_finish", mock_finish)
        resp = core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))

        assert core._store is None
        assert core._transport is None
        mock_heartbeat.stop.assert_called_once()
        mock_finish.assert_called_once()
        assert resp.success is False
        assert "saved locally" in resp.message


# ============================================================
# TestCorePythonGuard
# ============================================================


class TestCorePythonGuard:
    def test_double_start_raises(self, tmp_path):
        core = CorePython("local")
        core._ctx = make_core_ctx(tmp_path)
        core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        with pytest.raises(RuntimeError, match="already started"):
            core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

    def test_finish_before_start_raises(self, tmp_path):
        core = CorePython("local")
        core._ctx = make_core_ctx(tmp_path)

        with pytest.raises(RuntimeError, match="not started"):
            core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))

    def test_fork_raises(self, tmp_path):
        core = CorePython("disabled")

        with pytest.raises(RuntimeError, match="should not be called"):
            core.fork()
