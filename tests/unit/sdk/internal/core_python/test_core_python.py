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

from unittest.mock import MagicMock

import pytest

from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord, StartResponse
from swanlab.sdk.internal.context import RunConfig, RunContext
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.settings import Settings


def make_ctx(tmp_path, mode: str) -> RunContext:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    settings = Settings.model_validate({"mode": mode, "run": {"id": "test-run-id"}})
    return RunContext(RunConfig(run_dir=run_dir, settings=settings))


def make_start_record() -> StartRecord:
    return StartRecord(id="test-run-id", project="test-project")


# ============================================================
# TestCorePythonStart
# ============================================================


class TestCorePythonStart:
    @pytest.mark.parametrize("mode", ["disabled", "local", "offline"])
    def test_non_cloud_modes(self, tmp_path, mode):
        ctx = make_ctx(tmp_path, mode)
        core = CorePython(ctx)
        resp = core.deliver_run_start(make_start_record())

        assert resp.success is True
        assert core._transport is None
        if mode == "disabled":
            assert core._store is None
        else:
            assert core._store is not None

    def test_cloud_mode(self, tmp_path, monkeypatch):
        ctx = make_ctx(tmp_path, "cloud")
        core = CorePython(ctx)
        record = make_start_record()
        mock_deliver = MagicMock(return_value=StartResponse(success=True, message="OK", run=record))
        monkeypatch.setattr(core, "_report_run_start", mock_deliver)

        resp = core.deliver_run_start(record)

        assert resp.success is True
        assert core._store is not None
        assert core._transport is not None
        mock_deliver.assert_called_once_with(record)


# ============================================================
# TestCorePythonPublish
# ============================================================


class TestCorePythonPublish:
    def test_disabled_skips_silently(self, tmp_path):
        ctx = make_ctx(tmp_path, "disabled")
        core = CorePython(ctx)
        core.deliver_run_start(make_start_record())
        core.publish([])

    def test_not_started_skips_silently(self, tmp_path):
        ctx = make_ctx(tmp_path, "local")
        core = CorePython(ctx)
        core.publish([])


# ============================================================
# TestCorePythonFinish
# ============================================================


class TestCorePythonFinish:
    @pytest.mark.parametrize("mode", ["disabled", "local", "offline"])
    def test_non_cloud_modes(self, tmp_path, mode):
        ctx = make_ctx(tmp_path, mode)
        core = CorePython(ctx)
        core.deliver_run_start(make_start_record())

        resp = core.deliver_run_finish(FinishRecord())

        assert resp.success is True
        assert core._store is None

    def test_cloud_closes_store_and_transport(self, tmp_path, monkeypatch):
        ctx = make_ctx(tmp_path, "cloud")
        core = CorePython(ctx)

        mock_start = MagicMock(return_value=StartResponse(success=True, message="OK", run=make_start_record()))
        monkeypatch.setattr(core, "_report_run_start", mock_start)
        core.deliver_run_start(make_start_record())

        mock_finish = MagicMock(return_value=None)
        monkeypatch.setattr(core, "_report_run_finish", mock_finish)
        resp = core.deliver_run_finish(FinishRecord())

        assert core._store is None
        assert core._transport is None
        mock_finish.assert_called_once()
        assert resp.success is False
        assert "saved locally" in resp.message


# ============================================================
# TestCorePythonGuard
# ============================================================


class TestCorePythonGuard:
    def test_double_start_raises(self, tmp_path):
        ctx = make_ctx(tmp_path, "local")
        core = CorePython(ctx)
        core.deliver_run_start(make_start_record())

        with pytest.raises(RuntimeError, match="already started"):
            core.deliver_run_start(make_start_record())

    def test_finish_before_start_raises(self, tmp_path):
        ctx = make_ctx(tmp_path, "local")
        core = CorePython(ctx)

        with pytest.raises(RuntimeError, match="not started"):
            core.deliver_run_finish(FinishRecord())

    def test_fork_raises(self, tmp_path):
        ctx = make_ctx(tmp_path, "disabled")
        core = CorePython(ctx)

        with pytest.raises(RuntimeError, match="should not be called"):
            core.fork()
