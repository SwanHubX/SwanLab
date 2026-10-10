"""
@author: cunyue
@file: test_sync_helper.py.py
@time: 2026/5/17 19:30
@description: 测试sync的工具函数
"""

import pytest
from pydantic import ValidationError

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState, OperationStats
from swanlab.sdk.cmd.sync import ensure_run_dir, sync
from swanlab.sdk.internal.settings import Settings


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
    """Sync 在 CoreSyncPython.deliver_sync_flush 中校验凭证并自建 client，无 api_key 时阻断"""

    class FakeCoreSync:
        def deliver_sync_start(self, req):
            return type("Response", (), {"success": True, "message": "OK"})()

        def deliver_sync_flush(self):
            # 模拟 CoreSyncPython 在 deliver_sync_flush 中校验凭证失败的响应
            return type(
                "Response",
                (),
                {
                    "success": False,
                    "message": "Sync requires a valid API key. Login via `swanlab.login()` or configure SWANLAB_API_KEY / .netrc credentials.",
                },
            )()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: FakeCoreSync())

    with pytest.raises(RuntimeError, match="Sync requires a valid API key"):
        sync(tmp_path, settings=Settings())


def test_sync_logs_in_when_api_key_provided(tmp_path, monkeypatch):
    """Client 由 CoreSyncPython 自建与销毁（deliver_sync_flush 中创建，confirm_sync_finish 的 finally 中释放）"""

    class FakeCoreSync:
        def __init__(self):
            self.client_created = False
            self.client_released = False

        def deliver_sync_start(self, req):
            return DeliverSyncStartResponse(success=True, message="OK")

        def deliver_sync_flush(self):
            # 模拟 CoreSyncPython 在 deliver_sync_flush 中自建 client
            self.client_created = True
            return DeliverSyncFlushResponse(success=True, message="success", path="/username/project/run123")

        def confirm_sync_finish(self):
            # 模拟 CoreSyncPython 在 confirm_sync_finish 的 finally 中释放 client
            self.client_released = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)

    sync(tmp_path, settings=Settings(api_key="test-api-key", api_host="https://api.example.com"))

    assert core.client_created is True  # client 已在 deliver_sync_flush 中创建
    assert core.client_released is True  # client 已在 confirm_sync_finish 的 finally 中释放


def test_sync_rejected_while_run_active(tmp_path, monkeypatch):
    """进程内互斥：active run 期间 sync 在创建任何资源前被拒绝"""
    from swanlab.sdk.cmd.init import init
    from swanlab.sdk.internal.run import clear_run, has_run

    init(mode="disabled")
    assert has_run()
    try:
        with pytest.raises(RuntimeError, match="`swanlab.sync` requires no active Run"):
            sync(tmp_path, settings=Settings())
    finally:
        clear_run()


def test_sync_waits_for_progress_before_confirming(tmp_path, monkeypatch):
    """Sync 应当在进度展示后才调用 confirm_sync_finish"""

    class FakeCoreSync:
        def __init__(self):
            self.confirmed = False

        def deliver_sync_start(self, req):
            return DeliverSyncStartResponse(success=True, message="OK")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=True, message="success", path="/username/project/run123")

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

        def confirm_sync_finish(self):
            self.confirmed = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

    core = FakeCoreSync()

    def assert_not_confirmed_then_confirm(stats_fn, blocking_fn, **kwargs):
        # 进度展示期间 confirm_sync_finish 尚未被调用
        assert core.confirmed is False
        return blocking_fn()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr("swanlab.sdk.cmd.sync.run_with_progress", assert_not_confirmed_then_confirm)

    sync(tmp_path, settings=Settings(api_key="test-key"))

    assert core.confirmed is True  # 进度展示完成后 confirm_sync_finish 已被调用


def test_sync_uses_protocol_confirm_directly(tmp_path, monkeypatch):
    """Sync 直接调用 protocol 的 confirm_sync_finish 方法"""

    class FakeCoreSync:
        def __init__(self):
            self.confirmed = False

        def deliver_sync_start(self, req):
            return DeliverSyncStartResponse(success=True, message="OK")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=True, message="success", path="/username/project/run123")

        def get_operation_stats(self):
            return GetOperationStatsResponse(
                success=True, message="OK", stats=OperationStats(state=CoreState.CORE_STATE_FINISHED)
            )

        def confirm_sync_finish(self):
            self.confirmed = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

    core = FakeCoreSync()

    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)
    monkeypatch.setattr(
        "swanlab.sdk.cmd.sync.run_with_progress",
        lambda stats_fn, blocking_fn, **kwargs: blocking_fn(),
    )

    sync(tmp_path, settings=Settings(api_key="test-key"))
    assert core.confirmed is True


def test_sync_raises_when_deliver_sync_start_fails(tmp_path, monkeypatch):
    """deliver_sync_start 失败时应当抛出异常"""

    class FakeCoreSync:
        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=False, message="start failed")

        def deliver_sync_flush(self):
            raise AssertionError("flush should not be called")

    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: FakeCoreSync())

    with pytest.raises(RuntimeError, match="start failed"):
        sync(tmp_path, settings=Settings(api_key="test-key"))


def test_sync_raises_when_deliver_sync_flush_fails(tmp_path, monkeypatch):
    """deliver_sync_flush 失败时应当抛出异常，confirm_sync_finish 不应被调用"""

    class FakeCoreSync:
        def __init__(self):
            self.confirmed = False

        def deliver_sync_start(self, request):
            return DeliverSyncStartResponse(success=True, message="success")

        def deliver_sync_flush(self):
            return DeliverSyncFlushResponse(success=False, message="flush failed")

        def confirm_sync_finish(self):
            self.confirmed = True
            return ConfirmSyncFinishResponse(success=True, message="OK")

    core = FakeCoreSync()
    monkeypatch.setattr("swanlab.sdk.cmd.sync.impl.create_core_sync", lambda: core)

    with pytest.raises(RuntimeError, match="flush failed"):
        sync(tmp_path, settings=Settings(api_key="test-key"))

    assert core.confirmed is False  # flush 失败时 confirm 不应被调用
