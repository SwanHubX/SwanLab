"""
@author: cunyue
@file: test_sync_helper.py.py
@time: 2026/5/17 19:30
@description: 测试sync的工具函数
"""

import pytest
from pydantic import ValidationError

from swanlab.exceptions import AuthenticationError
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState, OperationStats
from swanlab.sdk.cmd.sync import build_sync_range, ensure_run_dir, format_sync_range, sync
from swanlab.sdk.internal.settings import Settings


def test_build_sync_range_returns_none_without_bounds():
    assert build_sync_range() is None


def test_build_sync_range_keeps_open_bounds():
    sync_range = build_sync_range(start_step=10, end_time="2026-10-20 12:00:00")

    assert sync_range is not None
    assert sync_range.HasField("start_step") is True
    assert sync_range.HasField("end_step") is False
    assert sync_range.HasField("start_time") is False
    assert sync_range.HasField("end_time") is True


def test_build_sync_range_formats_bounds():
    sync_range = build_sync_range(start_step=10, start_time="2026-10-20 10:00:00")

    assert sync_range is not None
    assert format_sync_range(sync_range) == "step 10..*, time 2026-10-20 10:00:00..*"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"start_step": -1},
        {"start_step": 10, "end_step": 5},
        {"start_time": "2026-10-20 12:00:00", "end_time": "2026-10-20 10:00:00"},
        {"start_time": "not-a-time"},
    ],
)
def test_build_sync_range_rejects_invalid_bounds(kwargs):
    with pytest.raises(ValueError):
        build_sync_range(**kwargs)


def test_sync_rejects_invalid_range_before_authentication(tmp_path):
    with pytest.raises(ValueError, match="start_step"):
        sync(tmp_path, settings=Settings(), start_step=10, end_step=5)


def test_sync_warns_when_only_step_range_is_given(tmp_path, monkeypatch):
    """step 范围不覆盖终端日志，静默丢弃前应先提示用户。"""
    warnings = []
    monkeypatch.setattr("swanlab.sdk.cmd.sync.console.warning", lambda message, *a, **k: warnings.append(message))
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: False)

    with pytest.raises(AuthenticationError):
        sync(tmp_path, settings=Settings(), start_step=0)

    assert len(warnings) == 1
    assert "log" in warnings[0].lower()


def test_ensure_run_dir_ok(tmp_path):
    assert ensure_run_dir(tmp_path) == tmp_path.resolve()


def test_ensure_run_dir_requires_directory(tmp_path):
    file_path = tmp_path / "run.txt"
    file_path.write_text("", encoding="utf-8")

    with pytest.raises(ValidationError):
        ensure_run_dir(file_path)


def test_ensure_run_dir_requires_readable(tmp_path, monkeypatch):
    monkeypatch.setattr("swanlab.sdk.cmd.sync.os.access", lambda *_: False)

    with pytest.raises(PermissionError):
        ensure_run_dir(tmp_path)


def test_sync_requires_api_key_when_client_missing(tmp_path, monkeypatch):
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: False)

    with pytest.raises(AuthenticationError, match="Please login first"):
        sync(tmp_path, settings=Settings())


def test_sync_logs_in_when_api_key_provided(tmp_path, monkeypatch):
    logged_in = False
    login_calls = []

    def client_exists():
        return logged_in

    def login_raw(api_key, *, host, print_welcome):
        nonlocal logged_in
        login_calls.append((api_key, host, print_welcome))
        logged_in = True

    class FakeCoreSync:
        def confirm_sync_finish(self):
            return type("Response", (), {"success": True, "message": "OK"})()

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", client_exists)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.login_raw", login_raw)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: FakeCoreSync())
    monkeypatch.setattr("swanlab.sdk.cmd.sync._deliver_sync_start", lambda *args, **kwargs: None)

    sync(tmp_path, settings=Settings(api_key="test-api-key", api_host="https://api.example.com"))

    assert login_calls == [("test-api-key", "https://api.example.com", False)]


def test_sync_waits_for_progress_before_confirming(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

        def confirm_sync_finish(self):
            self.confirmed = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

    core = FakeCoreSync()

    def assert_not_confirmed_then_confirm(stats_fn, blocking_fn, **kwargs):
        assert core.confirmed is False
        return blocking_fn()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._deliver_sync_start", lambda *args, **kwargs: None)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.run_with_progress", assert_not_confirmed_then_confirm)

    sync(tmp_path, settings=Settings())
    assert core.confirmed is True


def test_sync_uses_protocol_confirm_directly(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

        def confirm_sync_finish(self):
            self.confirmed = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

    core = FakeCoreSync()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._deliver_sync_start", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        "swanlab.sdk.cmd.sync.run_with_progress",
        lambda stats_fn, blocking_fn, **kwargs: blocking_fn(),
    )

    sync(tmp_path, settings=Settings())
    assert core.confirmed is True


def test_sync_raises_when_deliver_sync_start_fails(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=False, message="start failed")

        def deliver_sync_flush(self):
            raise AssertionError("flush should not be called")

        def confirm_sync_finish(self):
            self.confirmed = True
            return type("Response", (), {"success": True, "message": "OK"})()

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="start failed"):
        sync(tmp_path, settings=Settings())

    assert core.confirmed is False


def test_sync_raises_when_deliver_sync_flush_fails(tmp_path, monkeypatch):
    class FakeCoreSync:
        confirmed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=True, message="success")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=False, message="flush failed")

        def confirm_sync_finish(self):
            self.confirmed = True
            return type("Response", (), {"success": True, "message": "OK"})()

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="flush failed"):
        sync(tmp_path, settings=Settings())

    assert core.confirmed is False


def test_sync_fills_api_key_from_active_client_and_closes(tmp_path, monkeypatch):
    captured_settings = []

    class RecordingCoreSync:
        closed = False
        aborted = False

        def deliver_sync_start(self, request):
            captured_settings.append(request.core_settings)
            return DeliverSyncStartResponse(success=True, message="success")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=True, message="success", path="/@u/p/s")

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

        def confirm_sync_finish(self):
            return ConfirmSyncFinishResponse(success=True, message="OK")

        def close(self):
            self.closed = True

        def abort(self, cause=None):
            self.aborted = True

    core = RecordingCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.current_api_key", lambda: "runtime-secret-key")
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.cmd.sync.run_with_progress",
        lambda stats_fn, blocking_fn, **kwargs: blocking_fn(),
    )

    sync(tmp_path, settings=Settings())

    assert len(captured_settings) == 1
    assert captured_settings[0].api_key == "runtime-secret-key"
    assert core.closed is True
    assert core.aborted is False


def test_sync_aborts_core_on_exception_and_keyboard_interrupt(tmp_path, monkeypatch):
    for exc_to_raise in (RuntimeError("boom"), KeyboardInterrupt()):

        class FailingCoreSync:
            def __init__(self) -> None:
                self.closed = False
                self.aborted_causes = []

            def deliver_sync_start(self, request):
                raise exc_to_raise

            def deliver_sync_flush(self):
                return DeliverSyncFlushResponse(success=True, message="success")

            def get_operation_stats(self):
                return GetOperationStatsResponse(success=True)

            def confirm_sync_finish(self):
                return ConfirmSyncFinishResponse(success=True)

            def close(self):
                self.closed = True

            def abort(self, cause=None):
                self.aborted_causes.append(cause)

        core = FailingCoreSync()
        monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.sync.client.current_api_key", lambda: None)
        monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)

        with pytest.raises(type(exc_to_raise)):
            sync(tmp_path, settings=Settings())

        assert len(core.aborted_causes) == 1
        assert core.aborted_causes[0] is exc_to_raise
        assert core.closed is False


def test_sync_stops_polling_when_core_dies_mid_upload(tmp_path, monkeypatch):
    """core 在上传中退出：连续 stats 失败必须收敛，最终把异常抛给调用方并 abort 回收。"""

    class DyingCoreSync:
        def __init__(self) -> None:
            self.stats_calls = 0
            self.aborted_causes = []
            self.closed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=True, message="success")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=True, message="success")

        def get_operation_stats(self):
            self.stats_calls += 1
            if self.stats_calls == 1:
                return GetOperationStatsResponse(
                    success=True,
                    message="OK",
                    stats=OperationStats(
                        state=CoreState.CORE_STATE_RUNNING,
                        total_records=10,
                        uploaded_records=0,
                    ),
                )
            return GetOperationStatsResponse(success=False, message="[RPC/GRPC_UNAVAILABLE] connect failed")

        def confirm_sync_finish(self):
            raise RuntimeError("core is gone")

        def close(self):
            self.closed = True

        def abort(self, cause=None):
            self.aborted_causes.append(cause)

    core = DyingCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.exists", lambda: True)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.client.current_api_key", lambda: None)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr("swanlab.sdk.cmd.sync._deliver_sync_start", lambda *args, **kwargs: None)
    monkeypatch.setattr("swanlab.sdk.internal.run.progress._POLL_INTERVAL", 0)
    monkeypatch.setattr("swanlab.sdk.internal.run.progress._MAX_CONSECUTIVE_STATS_FAILURES", 2)

    with pytest.raises(RuntimeError, match="core is gone"):
        sync(tmp_path, settings=Settings())

    assert core.stats_calls <= 5
    assert len(core.aborted_causes) == 1
    assert core.closed is False
