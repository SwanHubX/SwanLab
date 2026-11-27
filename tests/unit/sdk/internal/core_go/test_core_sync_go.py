"""
@author: Nexisato
@file: test_core_sync_go.py
@time: 2026/11/18
@description: CoreSyncGo Owner Mode 生命周期、失败回滚与幂等关闭测试。
"""

from pathlib import Path
from typing import List, Optional, Tuple

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
    DeliverSyncStartResponse,
)
from swanlab.proto.swanlab.operation.v1.operation_pb2 import OperationStats
from swanlab.sdk.internal.core_go import CoreGoError, CoreSyncGo, ReasonCode, Stage
from swanlab.sdk.internal.core_go import sync as sync_module
from swanlab.sdk.internal.core_go.config import TEARDOWN_TIMEOUT
from swanlab.sdk.internal.core_go.sync import SyncState
from swanlab.sdk.protocol import is_terminal_stats_error


class FakeCoreRpc:
    instances: List["FakeCoreRpc"] = []

    def __init__(self, channel: object) -> None:
        self.channel = channel
        self.spinups = 0
        self.teardowns: List[float] = []
        self.spinup_error: Optional[BaseException] = None
        FakeCoreRpc.instances.append(self)

    def spinup(self) -> None:
        self.spinups += 1
        if self.spinup_error is not None:
            raise self.spinup_error

    def teardown(self, timeout: float = 0.0) -> None:
        self.teardowns.append(timeout)


class FakeSyncRpc:
    instances: List["FakeSyncRpc"] = []

    def __init__(self, channel: object) -> None:
        self.channel = channel
        self.start_calls: List[DeliverSyncStartRequest] = []
        self.flush_calls = 0
        self.stats_calls = 0
        self.confirm_calls = 0
        self.start_error: Optional[BaseException] = None
        self.flush_error: Optional[BaseException] = None
        self.stats_error: Optional[BaseException] = None
        self.confirm_error: Optional[BaseException] = None
        self.start_response = DeliverSyncStartResponse(success=True, message="success")
        self.flush_response = DeliverSyncFlushResponse(success=True, message="success", path="/@u/p/s")
        self.stats_response = GetOperationStatsResponse(
            success=True, stats=OperationStats(total_records=3, uploaded_records=3)
        )
        self.confirm_response = ConfirmSyncFinishResponse(success=True, message="success")
        FakeSyncRpc.instances.append(self)

    def deliver_sync_start(self, request: DeliverSyncStartRequest) -> DeliverSyncStartResponse:
        self.start_calls.append(request)
        if self.start_error is not None:
            raise self.start_error
        return self.start_response

    def deliver_sync_flush(self) -> DeliverSyncFlushResponse:
        self.flush_calls += 1
        if self.flush_error is not None:
            raise self.flush_error
        return self.flush_response

    def get_operation_stats(self) -> GetOperationStatsResponse:
        self.stats_calls += 1
        if self.stats_error is not None:
            raise self.stats_error
        return self.stats_response

    def confirm_sync_finish(self) -> ConfirmSyncFinishResponse:
        self.confirm_calls += 1
        if self.confirm_error is not None:
            raise self.confirm_error
        return self.confirm_response


@pytest.fixture
def patch_sync_dependencies(
    install_fake_launcher, monkeypatch
) -> Tuple[List[object], List[FakeCoreRpc], List[FakeSyncRpc]]:
    FakeCoreRpc.instances.clear()
    FakeSyncRpc.instances.clear()
    launchers, _ = install_fake_launcher(sync_module, exit_code=0)
    monkeypatch.setattr(sync_module, "CoreServiceRpc", FakeCoreRpc)
    monkeypatch.setattr(sync_module, "CoreSyncRpc", FakeSyncRpc)
    return launchers, FakeCoreRpc.instances, FakeSyncRpc.instances


def test_deliver_sync_start_without_binary_fails_and_closes():
    core = CoreSyncGo()

    resp = core.deliver_sync_start(DeliverSyncStartRequest())

    assert resp.success is False
    assert core.last_error is not None
    assert core.last_error.reason is ReasonCode.BINARY_NOT_FOUND
    assert core._state is SyncState.CLOSED  # noqa: SLF001


def test_full_owner_lifecycle_and_idempotent_close(patch_sync_dependencies, tmp_path: Path):
    launchers, core_rpcs, sync_rpcs = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    start_resp = core.deliver_sync_start(DeliverSyncStartRequest(project="proj"))
    assert start_resp.success is True
    assert launchers[0].detach is False
    assert launchers[0].started is True
    assert core_rpcs[0].spinups == 1
    assert core._state is SyncState.STARTED  # noqa: SLF001

    flush_resp = core.deliver_sync_flush()
    assert flush_resp.success is True
    assert flush_resp.path == "/@u/p/s"
    assert core._state is SyncState.RUNNING  # noqa: SLF001

    stats_resp = core.get_operation_stats()
    assert stats_resp.success is True
    assert stats_resp.stats.total_records == 3

    confirm_resp = core.confirm_sync_finish()
    assert confirm_resp.success is True
    assert core._state is SyncState.FINISHED  # noqa: SLF001

    core.close()
    core.close()
    core.abort()

    assert core_rpcs[0].teardowns == [TEARDOWN_TIMEOUT]
    assert launchers[0].closed == 1
    assert launchers[0].aborted == 0
    assert core.exit_code == 0
    assert core._state is SyncState.CLOSED  # noqa: SLF001
    assert sync_rpcs[0].confirm_calls == 1


def test_deliver_sync_start_rejection_aborts_owner_process(patch_sync_dependencies, monkeypatch, tmp_path: Path):
    launchers, core_rpcs, _ = patch_sync_dependencies

    class RejectingSyncRpc(FakeSyncRpc):
        def __init__(self, channel: object) -> None:
            super().__init__(channel)
            self.start_error = CoreGoError(Stage.RUN_START, ReasonCode.RUN_START_REJECTED, "bad run dir")

    FakeSyncRpc.instances.clear()
    monkeypatch.setattr(sync_module, "CoreSyncRpc", RejectingSyncRpc)
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    resp = core.deliver_sync_start(DeliverSyncStartRequest())

    assert resp.success is False
    assert "bad run dir" in resp.message
    assert core.last_error is not None
    assert core.last_error.reason is ReasonCode.RUN_START_REJECTED
    assert core_rpcs[0].teardowns == [TEARDOWN_TIMEOUT]
    assert launchers[0].aborted == 1
    assert core._state is SyncState.CLOSED  # noqa: SLF001


def test_deliver_sync_flush_failure_aborts_and_returns_failed_response(patch_sync_dependencies, tmp_path: Path):
    launchers, core_rpcs, sync_rpcs = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    # 未 start 即 flush，不拉起进程，直接拒绝
    pre_resp = core.deliver_sync_flush()
    assert pre_resp.success is False
    assert core.last_error is not None
    assert core.last_error.reason is ReasonCode.NOT_STARTED

    assert core.deliver_sync_start(DeliverSyncStartRequest()).success is True
    sync_rpcs[0].flush_error = CoreGoError(Stage.RUN_START, ReasonCode.RUN_START_REJECTED, "401 unauthorized")

    flush_resp = core.deliver_sync_flush()

    assert flush_resp.success is False
    assert "401 unauthorized" in flush_resp.message
    assert core_rpcs[0].teardowns == [TEARDOWN_TIMEOUT]
    assert launchers[0].aborted == 1
    assert core._state is SyncState.CLOSED  # noqa: SLF001


def test_keyboard_interrupt_aborts_and_reraises(patch_sync_dependencies, tmp_path: Path):
    launchers, core_rpcs, sync_rpcs = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    assert core.deliver_sync_start(DeliverSyncStartRequest()).success is True
    assert core.deliver_sync_flush().success is True
    sync_rpcs[0].confirm_error = KeyboardInterrupt()

    with pytest.raises(KeyboardInterrupt):
        core.confirm_sync_finish()

    assert core_rpcs[0].teardowns == [TEARDOWN_TIMEOUT]
    assert launchers[0].aborted == 1
    assert core._state is SyncState.CLOSED  # noqa: SLF001


def test_get_operation_stats_handles_prestart_and_transient_errors(patch_sync_dependencies, tmp_path: Path):
    _, _, sync_rpcs = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    pre_stats = core.get_operation_stats()
    assert pre_stats.success is False
    assert is_terminal_stats_error(pre_stats) is True

    assert core.deliver_sync_start(DeliverSyncStartRequest()).success is True
    assert core.deliver_sync_flush().success is True

    sync_rpcs[0].stats_error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "temporary blip")
    err_stats = core.get_operation_stats()
    assert err_stats.success is False
    assert "temporary blip" in err_stats.message
    # 进度查询失败不应中止 RUNNING 会话
    assert core._state is SyncState.RUNNING  # noqa: SLF001

    # Confirm 过程中并发 stats 依然允许查询
    core._state = SyncState.FINISHING  # noqa: SLF001
    sync_rpcs[0].stats_error = None
    assert core.get_operation_stats().success is True


def test_confirm_sync_finish_rpc_error_returns_failed_response(patch_sync_dependencies, tmp_path: Path):
    launchers, core_rpcs, sync_rpcs = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")

    assert core.deliver_sync_start(DeliverSyncStartRequest()).success is True
    assert core.deliver_sync_flush().success is True
    sync_rpcs[0].confirm_error = CoreGoError(Stage.RPC, ReasonCode.RUN_FINISH_REJECTED, "Failed to finish run")

    resp = core.confirm_sync_finish()

    assert resp.success is False
    assert "Failed to finish run" in resp.message
    assert core._state is SyncState.FINISHED  # noqa: SLF001
    core.close()
    assert core_rpcs[0].teardowns == [TEARDOWN_TIMEOUT]
    assert launchers[0].closed == 1


def test_fork_child_does_not_teardown_parent_core(patch_sync_dependencies, monkeypatch, tmp_path: Path):
    launchers, core_rpcs, _ = patch_sync_dependencies
    core = CoreSyncGo(binary=tmp_path / "swanlab-core")
    assert core.deliver_sync_start(DeliverSyncStartRequest()).success is True

    from swanlab.sdk.internal.pkg import fork as pkg_fork

    monkeypatch.setattr(pkg_fork, "is_forked", lambda _owner_pid: True)
    core.close()
    core.abort()

    assert core_rpcs[0].teardowns == []
    assert launchers[0].closed == 0
    assert launchers[0].aborted == 0
