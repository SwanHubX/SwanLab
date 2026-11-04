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
  - TestCorePythonSkipStore: core.skip_store 主链路（无本地 datastore、Null watcher + live save 回归）
"""

from unittest.mock import MagicMock, patch

import pytest
import requests

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    CORE_ERROR_KIND_AUTH,
    CORE_ERROR_KIND_OTHER,
    CORE_ERROR_KIND_TRANSIENT,
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy, SaveRecord
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings as CoreSettingsPb
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.store import DataStoreWriter, NullDataStoreWriter
from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
from swanlab.sdk.internal.core_python.watcher import FileWatcher, NullFileWatcher


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

    @pytest.mark.parametrize("new_experiment", [True, False])
    def test_report_run_start_uses_startup_retries(self, tmp_path, monkeypatch, new_experiment):
        """online 启动握手应使用统一的 STARTUP_RETRIES。"""
        from swanlab.sdk.internal.core_python import core as core_module
        from swanlab.sdk.internal.core_python.client import STARTUP_RETRIES
        from swanlab.sdk.internal.core_python.utils import PrepareExperimentStartResult

        core = CorePython("online")
        core._ctx = make_core_ctx(tmp_path)
        request = make_start_request(tmp_path, make_start_record())
        prepare = MagicMock(
            return_value=PrepareExperimentStartResult(
                username="test-user",
                project="test-project",
                project_data={
                    "cuid": "test-project-id",
                    "name": "test-project",
                    "path": "/test-user/test-project",
                    "visibility": "PRIVATE",
                    "group": {"username": "test-user"},
                    "_count": {"experiments": 0, "contributors": 0, "collaborators": 0, "clones": 0},
                },
                experiment={
                    "cuid": "test-experiment-id",
                    "slug": "test-run-id",
                    "name": "test-run",
                    "createdAt": "2024-07-01T00:00:00Z",
                },
                new_experiment=new_experiment,
                name="test-run",
                color="#000000",
            )
        )
        summary = MagicMock(return_value=None)
        monkeypatch.setattr(core_module, "prepare_experiment_start", prepare)
        monkeypatch.setattr(core_module, "get_experiment_summary", summary)

        core._report_run_start(request.start_record)

        assert prepare.call_args.kwargs["retries"] == STARTUP_RETRIES
        if new_experiment:
            summary.assert_not_called()
        else:
            assert summary.call_args.kwargs["retries"] == STARTUP_RETRIES

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
# TestCorePythonOnlineStartFailure
# ============================================================


class TestCorePythonOnlineStartFailure:
    def test_connection_error_returns_transient(self, tmp_path, monkeypatch):
        core = CorePython("online")
        monkeypatch.setattr(core, "_report_run_start", MagicMock(side_effect=requests.ConnectionError("boom")))

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_TRANSIENT
        assert "boom" in resp.message
        assert core._active is False

    def test_authentication_error_returns_auth(self, tmp_path, monkeypatch):
        core = CorePython("online")
        monkeypatch.setattr(
            core,
            "_report_run_start",
            MagicMock(side_effect=AuthenticationError("invalid API key")),
        )

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_AUTH

    def test_api_503_returns_transient(self, tmp_path, monkeypatch):
        core = CorePython("online")
        error = ApiError(
            MagicMock(status_code=503, url="https://api.test/project"),
            method="GET",
            code=503,
            message="unavailable",
        )
        monkeypatch.setattr(core, "_report_run_start", MagicMock(side_effect=error))

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_TRANSIENT

    def test_business_error_returns_other(self, tmp_path, monkeypatch):
        core = CorePython("online")
        monkeypatch.setattr(core, "_report_run_start", MagicMock(side_effect=RuntimeError("experiment not found")))

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_OTHER

    def test_malformed_response_returns_transient(self, tmp_path, monkeypatch):
        """200 响应体畸形（decode_response 退化为 text）表现为 TypeError，允许降级。"""
        core = CorePython("online")
        monkeypatch.setattr(
            core, "_report_run_start", MagicMock(side_effect=TypeError("string indices must be integers"))
        )

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_TRANSIENT

    @pytest.mark.parametrize(
        "error",
        [
            ValueError("Backend did not return createdAt; please upgrade swanlab-server"),
            AssertionError("experiment without name"),
        ],
    )
    def test_config_or_internal_error_returns_other(self, tmp_path, monkeypatch, error):
        """配置 / 后端兼容性 / 内部断言等确定性错误不允许降级。"""
        core = CorePython("online")
        monkeypatch.setattr(core, "_report_run_start", MagicMock(side_effect=error))

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_OTHER

    def test_local_start_failure_is_other_and_cleans_resources(self, tmp_path, monkeypatch):
        core = CorePython("online")

        def _report_run_start_and_set_attrs(record):
            set_online_params(core)
            core._metrics = MagicMock()
            return DeliverRunStartResponse(success=True, message="OK", run=make_start_record())

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)
        monkeypatch.setattr(
            "swanlab.sdk.internal.core_python.core.Transport",
            MagicMock(side_effect=RuntimeError("local boom")),
        )

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        assert resp.error_kind == CORE_ERROR_KIND_OTHER
        assert core._store is None
        assert core._transport is None
        assert core._active is False

    def test_start_failure_does_not_enqueue_transport_records(self, tmp_path, monkeypatch):
        """启动阶段只允许写本地 start record；失败清理时 transport 不应收到任何记录。"""
        core = CorePython("online")

        def _report_run_start_and_set_attrs(record):
            set_online_params(core)
            core._metrics = MagicMock()
            return DeliverRunStartResponse(success=True, message="OK", run=make_start_record())

        transport = MagicMock()
        transport.pending_records.return_value = []
        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.Transport", MagicMock(return_value=transport))
        monkeypatch.setattr(
            "swanlab.sdk.internal.core_python.core.Heartbeat",
            MagicMock(side_effect=RuntimeError("heartbeat boom")),
        )

        resp = core.deliver_run_start(make_start_request(tmp_path, make_start_record()))

        assert resp.success is False
        transport.put.assert_not_called()
        transport.request_finish.assert_called_once()
        transport.join.assert_called_once_with(timeout=1.0)
        assert core._transport is None
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
        monkeypatch.setattr("swanlab.sdk.internal.core_python.transport.sender.upload_log", lambda *a, **kw: None)
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

    @pytest.mark.parametrize(
        "summary",
        [
            "Upload did not finish: 2 record(s) retained (worker exited unexpectedly).",
            "Upload did not finish: 1 record(s) retained (last failure: swanlab_api - http 503).",
        ],
        ids=["worker-exited", "last-failure"],
    )
    def test_confirm_run_finish_failure_keeps_transport_and_stops_components(self, tmp_path, monkeypatch, summary):
        """上传未排空时保留诊断数据，停止心跳并返回失败摘要。"""
        core = CorePython("online")
        core._ctx = make_core_ctx(tmp_path)
        set_online_params(core)
        transport = MagicMock()
        transport.finish.return_value = False
        transport.finish_summary.return_value = summary
        heartbeat = MagicMock()
        tracker = UploadTracker()
        tracker.set_state(CoreState.CORE_STATE_RUNNING)
        stop_experiment_mock = MagicMock()
        core._transport = transport
        core._heartbeat = heartbeat
        core._tracker = tracker
        core._pending_online_finish_record = FinishRecord()
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", stop_experiment_mock)

        resp = core.confirm_run_finish()

        assert isinstance(resp, ConfirmRunFinishResponse)
        assert resp.success is False
        assert resp.message == summary
        assert core._transport is transport
        transport.finish.assert_called_once_with(timeout=None)
        assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED
        stop_experiment_mock.assert_not_called()

    def test_confirm_run_finish_returns_failure_when_report_fails(self, tmp_path, monkeypatch):
        """后端上报失败时 confirm 释放传输组件并返回失败，本地数据已保存。"""
        core = CorePython("online")
        core._ctx = make_core_ctx(tmp_path)
        set_online_params(core)
        tracker = UploadTracker()
        tracker.set_state(CoreState.CORE_STATE_RUNNING)
        transport = MagicMock()
        transport.finish.return_value = True
        heartbeat = MagicMock()
        core._tracker = tracker
        core._transport = transport
        core._heartbeat = heartbeat
        core._pending_online_finish_record = FinishRecord()

        def fail_report(*args, **kwargs):
            raise ConnectionError("network unreachable")

        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", fail_report)

        resp = core.confirm_run_finish()

        assert isinstance(resp, ConfirmRunFinishResponse)
        assert resp.success is False
        assert "Failed to finish run" in resp.message
        assert core._transport is None
        heartbeat.stop.assert_called_once()
        assert tracker.snapshot().state == CoreState.CORE_STATE_FINISHED


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


# ============================================================
# TestCorePythonCleanup
# ============================================================


class TestCorePythonCleanup:
    def test_transport_is_stopped_even_with_pending_records(self):
        """启动失败清理时即使缓冲非空，也必须通知 transport 退出，避免线程泄漏。"""
        core = CorePython("online")
        transport = MagicMock()
        transport.pending_records.return_value = [MagicMock()]
        core._transport = transport

        core._cleanup_failed_start()

        transport.request_finish.assert_called_once()
        transport.join.assert_called_once_with(timeout=1.0)
        assert core._transport is None

    def test_watcher_stops_when_store_close_fails(self):
        """store 关闭失败不能阻止 watcher 停止，且首个异常仍向调用方抛出。"""
        core = CorePython("offline")
        store = MagicMock()
        store.close.side_effect = OSError("disk error")
        watcher = MagicMock()
        core._store = store
        core._watcher = watcher

        with pytest.raises(OSError, match="disk error"):
            core._release_local_storage()

        store.close.assert_called_once()
        watcher.stop.assert_called_once()
        assert core._store is None
        assert core._watcher is None


# ============================================================
# TestCorePythonSkipStore
# ============================================================


class TestCorePythonSkipStore:
    """core.skip_store=true（仅 online 合法）：本地不产生 datastore，记录仍进 Transport。"""

    def _start_online_core(self, tmp_path, monkeypatch, skip_store: bool) -> CorePython:
        core = CorePython("online")
        record = make_start_record()
        mock_start = MagicMock(return_value=DeliverRunStartResponse(success=True, message="OK", run=record))

        def _report_run_start_and_set_attrs(rec):
            set_online_params(core)
            core._metrics = MagicMock()
            return mock_start(rec)

        monkeypatch.setattr(core, "_report_run_start", _report_run_start_and_set_attrs)
        req = make_start_request(tmp_path, record)
        req.core_settings.skip_store = skip_store
        with patch("swanlab.sdk.internal.core_python.core.Heartbeat"):
            resp = core.deliver_run_start(req)
        assert resp.success is True
        return core

    def test_skip_store_creates_no_datastore_file(self, tmp_path, monkeypatch):
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=True)
        # from_proto 已解析该字段，skip 设置后 store 不产生文件
        assert core._ctx.config.skip_store is True
        assert core._store is not None
        assert isinstance(core._store, NullDataStoreWriter)
        assert isinstance(core._watcher, NullFileWatcher)
        assert not core._ctx.run_file.exists()
        assert core._transport is not None

    def test_default_creates_datastore_file(self, tmp_path, monkeypatch):
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=False)
        assert core._ctx.config.skip_store is False
        assert core._store is not None
        assert isinstance(core._store, DataStoreWriter)
        assert isinstance(core._watcher, FileWatcher)
        assert core._ctx.run_file.exists()

    def test_live_save_still_reaches_transport(self, tmp_path, monkeypatch):
        """skip 设置时 _on_file_changed 不能早退，record 必须进 Transport。"""
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=True)
        transport = MagicMock()
        core._transport = transport

        core._on_file_changed(SaveRecord())

        transport.put.assert_called_once()
        assert not core._ctx.run_file.exists()

    def test_live_save_after_finish_is_ignored(self, tmp_path, monkeypatch):
        """finish 将 store 置 None 后，watcher 回调应被早退条件拦住，不触发 assert。"""
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=True)
        monkeypatch.setattr("swanlab.sdk.internal.core_python.core.stop_experiment", lambda *a, **kw: None)
        # _store_finish 在非 FINISHED 状态会生成 error log record 进 transport，需让上传成功以便排空
        monkeypatch.setattr("swanlab.sdk.internal.core_python.transport.sender.upload_log", lambda *a, **kw: None)
        core.deliver_run_finish(DeliverRunFinishRequest(finish_record=FinishRecord()))
        assert core._store is None

        core._on_file_changed(SaveRecord())  # 不抛 assert

        assert core._store is None
        assert not core._ctx.run_file.exists()

    def test_store_records_skips_serialization(self, tmp_path, monkeypatch):
        """skip 下 _store_records 直接返回，不调用 SerializeToString。"""
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=True)

        record = MagicMock()
        core._store_records([record, record, record])

        record.SerializeToString.assert_not_called()

    def test_custom_save_skips_local_links(self, tmp_path, monkeypatch):
        """skip 下不创建镜像软链接、不填 target_path，也不触碰 files 目录；live 记录无监听路径，正常上传。"""
        core = self._start_online_core(tmp_path, monkeypatch, skip_store=True)
        transport = MagicMock()
        core._transport = transport
        assert isinstance(core._watcher, NullFileWatcher)
        source = tmp_path / "checkpoints" / "model.pt"
        source.parent.mkdir()
        source.write_bytes(b"weights")

        save = SaveRecord(
            name="checkpoints/model.pt",
            source_path=str(source),
            policy=SavePolicy.SAVE_POLICY_LIVE,
        )
        core.upsert_saves([save])

        assert save.target_path == ""
        assert source.read_bytes() == b"weights"  # 源文件只读不动
        assert not (core._ctx.config.run_dir / "files").exists()
        transport.put.assert_called_once()
