"""
@author: Nexisato
@file: test_core.py
@time: 2026/11/9
@description: CoreGo upsert 错误传播与协议层语义测试。

关键契约：
- CoreGo 把 upsert 失败抛给 consumer，并带上 retryable 分类；
- 默认 CoreProtocol 实现保持“记录日志后吞掉”的历史语义；
- CorePython 通过基类走的路径不受影响。
"""

from typing import List, Optional, cast

import pytest

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    CORE_ERROR_KIND_AUTH,
    CORE_ERROR_KIND_OTHER,
    CORE_ERROR_KIND_TRANSIENT,
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord, ScalarRecord
from swanlab.sdk.internal.core_go import core as core_module
from swanlab.sdk.internal.core_go.config import TEARDOWN_TIMEOUT, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.core import CoreGo, CoreState, DetachedCoreGo
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc
from swanlab.sdk.internal.pkg import console as pkg_console
from swanlab.sdk.protocol import CoreProtocol, CoreProtocolError, NullCoreProtocol, is_terminal_stats_error


class FakeRpc:
    """只实现 upsert_scalars / upsert_histograms，按需抛出预置异常。"""

    def __init__(self, error: Optional[BaseException] = None) -> None:
        self.error = error
        self.calls: List[tuple] = []
        self.histogram_calls: List[tuple] = []
        self.teardown_calls: List[float] = []

    def upsert_scalars(self, scalars: List[ScalarRecord], session_id: str) -> None:
        self.calls.append((scalars, session_id))
        if self.error is not None:
            raise self.error

    def upsert_histograms(self, histograms: List[HistogramRecord], session_id: str) -> None:
        self.histogram_calls.append((histograms, session_id))
        if self.error is not None:
            raise self.error

    def teardown(self, timeout: float = 0) -> None:
        self.teardown_calls.append(timeout)


def _running_core(rpc: FakeRpc) -> CoreGo:
    core = CoreGo("online")
    core._rpc = cast(CoreServiceRpc, rpc)  # noqa: SLF001 - 测试直接注入 fake RPC
    core._session_id = "session-1"  # noqa: SLF001
    core._state = CoreState.RUNNING  # noqa: SLF001
    return core


def test_core_go_opts_in_to_error_propagation():
    assert CoreProtocol.raises_on_upsert_error is False
    assert CoreGo.raises_on_upsert_error is True


def test_retryable_error_propagates():
    error = CoreGoError(Stage.RPC, ReasonCode.RESOURCE_EXHAUSTED, "server busy")
    core = _running_core(FakeRpc(error))

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value is error
    assert caught.value.retryable is True
    assert core.last_error is error


def test_terminal_error_propagates():
    error = CoreGoError(Stage.RPC, ReasonCode.INVALID_ARGUMENT, "bad record")
    core = _running_core(FakeRpc(error))

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.retryable is False


def test_histogram_upsert_failure_propagates_like_other_channels():
    """直方图与其余五类一致：失败必须抛给 consumer，而不是被协议层吞掉。"""
    error = CoreGoError(Stage.RPC, ReasonCode.RESOURCE_EXHAUSTED, "batch too large")
    rpc = FakeRpc(error)
    core = _running_core(rpc)
    record = HistogramRecord(key="weights", step=1)

    with pytest.raises(CoreGoError) as caught:
        core.upsert_histograms([record])

    assert caught.value is error
    assert caught.value.retryable is True
    assert rpc.histogram_calls == [([record], "session-1")]
    assert core.last_error is error


def test_unknown_exception_is_wrapped_as_terminal():
    core = _running_core(FakeRpc(RuntimeError("socket closed")))

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.reason is ReasonCode.GRPC_ERROR
    assert caught.value.retryable is False
    assert isinstance(caught.value.__cause__, RuntimeError)


def test_upsert_before_start_is_terminal_and_skips_rpc():
    rpc = FakeRpc()
    core = CoreGo("online")
    core._rpc = cast(CoreServiceRpc, rpc)  # noqa: SLF001

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.reason is ReasonCode.NOT_STARTED
    assert caught.value.retryable is False
    assert rpc.calls == []


def test_upsert_after_close_is_terminal():
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._state = CoreState.CLOSED  # noqa: SLF001

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.reason is ReasonCode.CLOSED
    assert caught.value.retryable is False
    assert rpc.calls == []


def test_core_go_does_not_log_on_upsert_failure(monkeypatch):
    """告警由 consumer 统一处理，CoreGo 不再自己打印，避免重复日志。"""
    called = []
    monkeypatch.setattr(pkg_console, "warning", lambda *args, **kwargs: called.append(args))

    core = _running_core(FakeRpc(CoreGoError(Stage.RPC, ReasonCode.RESOURCE_EXHAUSTED, "busy")))
    with pytest.raises(CoreGoError):
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert called == []


def test_default_protocol_swallows_upsert_error(monkeypatch):
    """默认实现（CorePython 路径）保持吞异常语义。"""
    monkeypatch.setattr(pkg_console, "trace", lambda *args, **kwargs: None)

    class BoomCore(NullCoreProtocol):
        def _upsert_scalars_when_online(self, scalars: List[ScalarRecord]) -> None:
            raise RuntimeError("boom")

    BoomCore("online").upsert_scalars([ScalarRecord(key="loss", step=1)])


def test_default_protocol_swallows_histogram_error(monkeypatch):
    """默认实现（CorePython 路径）对直方图同样吞异常。"""
    monkeypatch.setattr(pkg_console, "trace", lambda *args, **kwargs: None)

    class BoomCore(NullCoreProtocol):
        def _upsert_histograms_when_online(self, histograms: List[HistogramRecord]) -> None:
            raise RuntimeError("boom")

    BoomCore("online").upsert_histograms([HistogramRecord(key="weights", step=1)])


def test_core_protocol_error_carries_retryable():
    assert CoreProtocolError("retry", retryable=True).retryable is True
    assert CoreProtocolError("stop", retryable=False).retryable is False


@pytest.mark.parametrize(
    "message, expected",
    [
        ("[RPC/SESSION_NOT_FOUND] gone (grpc=NOT_FOUND)", True),
        ("[RPC/NOT_STARTED] run is not running", True),
        ("[RPC/CLOSED] channel closed", True),
        ("[RPC/GRPC_UNAVAILABLE] connect failed", False),
        ("plain failure", False),
    ],
)
def test_is_terminal_stats_error(message, expected):
    response = GetOperationStatsResponse(success=False, message=message)

    assert is_terminal_stats_error(response) is expected


def test_successful_stats_is_not_terminal():
    assert is_terminal_stats_error(GetOperationStatsResponse(success=True)) is False


def test_close_before_start_is_idempotent():
    core = CoreGo("online")

    core.close()
    core.close()

    assert core._state is CoreState.CLOSED  # noqa: SLF001


def test_close_after_finish_tears_down_and_releases_rpc():
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._state = CoreState.FINISHED  # noqa: SLF001

    core.close()

    assert core._state is CoreState.CLOSED  # noqa: SLF001
    assert core._rpc is None  # noqa: SLF001
    assert rpc.teardown_calls == [TEARDOWN_TIMEOUT]


def test_fork_returns_detached_placeholder(tmp_path):
    binary = tmp_path / "swanlab-core"
    core = CoreGo("online", binary=binary)

    detached = core.fork()

    assert isinstance(detached, DetachedCoreGo)
    assert detached is not core
    assert detached.fork() is detached
    assert detached._binary == binary  # noqa: SLF001
    detached.close()


def test_detached_core_never_touches_parent_lifecycle():
    detached = CoreGo("online").fork()

    assert detached.deliver_run_start(DeliverRunStartRequest()).success is False
    assert detached.deliver_run_finish(DeliverRunFinishRequest()).success is False
    assert detached.confirm_run_finish().success is False
    assert detached.get_operation_stats().success is False
    detached.upsert_scalars([ScalarRecord(key="loss", step=1)])
    detached.close()
    detached.abort()

    error = cast(CoreGoError, detached.last_error)
    assert error.reason is ReasonCode.FORK_UNSUPPORTED
    assert error.stage is Stage.FORK


def test_start_modes_share_enabled_path(monkeypatch):
    sentinel = DeliverRunStartResponse(success=True, session_id="s")
    calls = []
    monkeypatch.setattr(CoreGo, "_start_when_enabled", lambda self, req: calls.append(self._mode) or sentinel)

    for mode in ("online", "local", "offline"):
        assert CoreGo(mode).deliver_run_start(DeliverRunStartRequest()) is sentinel
    assert calls == ["online", "local", "offline"]


def test_finish_modes_share_enabled_path(monkeypatch):
    sentinel = DeliverRunFinishResponse(success=True)
    monkeypatch.setattr(CoreGo, "_finish_when_enabled", lambda self, req: sentinel)

    for mode in ("online", "local", "offline"):
        assert CoreGo(mode).deliver_run_finish(DeliverRunFinishRequest()) is sentinel


def test_start_without_binary_is_terminal():
    core = CoreGo("local")

    response = core.deliver_run_start(DeliverRunStartRequest())

    assert response.success is False
    assert "BINARY_NOT_FOUND" in response.message
    assert response.error_kind == CORE_ERROR_KIND_OTHER


def test_disabled_start_is_noop_without_binary():
    core = CoreGo("disabled")

    response = core.deliver_run_start(DeliverRunStartRequest())

    assert response.success is True


@pytest.mark.parametrize(
    "reason",
    [ReasonCode.GRPC_UNAVAILABLE, ReasonCode.DEADLINE_EXCEEDED, ReasonCode.CHANNEL_TIMEOUT],
)
def test_owned_core_process_exit_makes_any_retryable_error_terminal(fake_launcher_cls, reason):
    """独占进程已退出时，所有可重试错误都止损，不必只看 UNAVAILABLE。"""
    core = _running_core(FakeRpc(CoreGoError(Stage.RPC, reason, "rpc failed")))
    launcher = fake_launcher_cls(owns_process=True)
    launcher.process_exited = True
    core._launcher = cast(CoreLauncher, launcher)  # noqa: SLF001

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.reason is reason
    assert caught.value.retryable is False
    assert "process exited" in caught.value.message


def test_owned_core_process_exit_keeps_terminal_error_semantics(fake_launcher_cls):
    """进程退出只收窄可重试错误，非重试错误保持原始语义。"""
    error = CoreGoError(Stage.RPC, ReasonCode.INVALID_ARGUMENT, "bad record")
    core = _running_core(FakeRpc(error))
    launcher = fake_launcher_cls(owns_process=True)
    launcher.process_exited = True
    core._launcher = cast(CoreLauncher, launcher)  # noqa: SLF001

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value is error
    assert caught.value.retryable is False
    assert "process exited" not in caught.value.message


@pytest.mark.parametrize("action, initial_state", [("close", CoreState.FINISHED), ("abort", CoreState.RUNNING)])
def test_detached_mode_close_and_abort_skip_teardown_and_are_idempotent(fake_launcher_cls, action, initial_state):
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._detach = True  # noqa: SLF001
    launcher = fake_launcher_cls()
    core._launcher = cast(CoreLauncher, launcher)  # noqa: SLF001
    core._state = initial_state  # noqa: SLF001

    if action == "close":
        core.close()
        core.close()
        assert launcher.closed == 1
        assert launcher.aborted == 0
    else:
        core.abort(RuntimeError("boom"))
        core.abort()
        assert launcher.aborted == 1
        assert launcher.closed == 0

    assert rpc.teardown_calls == []
    assert core.exit_code is None
    assert core._state is CoreState.CLOSED  # noqa: SLF001


def test_owner_mode_abort_sends_teardown():
    rpc = FakeRpc()
    core = _running_core(rpc)

    core.abort(RuntimeError("boom"))

    assert rpc.teardown_calls == [TEARDOWN_TIMEOUT]
    assert core._state is CoreState.CLOSED  # noqa: SLF001


def test_unowned_launcher_skips_teardown_on_close_even_when_detach_false(fake_launcher_cls):
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._detach = False  # noqa: SLF001
    launcher = fake_launcher_cls(owns_process=False)
    core._launcher = cast(CoreLauncher, launcher)  # noqa: SLF001

    core.close()

    assert rpc.teardown_calls == []
    assert launcher.closed == 1
    assert core._state is CoreState.CLOSED  # noqa: SLF001


@pytest.mark.parametrize("first_mode", ["detached", "owner"])
def test_mixed_owner_and_detached_close_does_not_affect_peer(install_fake_launcher, monkeypatch, tmp_path, first_mode):
    binary = tmp_path / "swanlab-core"
    rpcs = {}

    launchers, _ = install_fake_launcher(
        core_module,
        exit_code=lambda detach: None if detach else 0,
        channel_factory=lambda _index, detach: "detached" if detach else "owner",
        abort_increments_closed=True,
    )

    class TaggedRpc(FakeRpc):
        def __init__(self, channel) -> None:
            super().__init__()
            self.channel = channel
            self.session_id = f"sid-{channel}"
            rpcs[channel] = self

        def spinup(self) -> None:
            return None

        def deliver_run_start(self, _req) -> DeliverRunStartResponse:
            return DeliverRunStartResponse(success=True, session_id=self.session_id)

    monkeypatch.setattr(core_module, "CoreServiceRpc", TaggedRpc)

    detached_core = CoreGo("online", binary=binary, detach=True, idle_timeout=2.0)
    owner_core = CoreGo("online", binary=binary, detach=False)

    if first_mode == "detached":
        assert detached_core.deliver_run_start(DeliverRunStartRequest()).success is True
        assert owner_core.deliver_run_start(DeliverRunStartRequest()).success is True
    else:
        assert owner_core.deliver_run_start(DeliverRunStartRequest()).success is True
        assert detached_core.deliver_run_start(DeliverRunStartRequest()).success is True

    launcher_kwargs = {
        launcher.detach: {"binary": launcher.binary, "idle_timeout": launcher.idle_timeout} for launcher in launchers
    }
    assert launcher_kwargs == {
        True: {"binary": binary, "idle_timeout": 2.0},
        False: {"binary": binary, "idle_timeout": core_module.DETACH_IDLE_TIMEOUT},
    }

    # Owner 退出只 teardown 自己的 Core，不得影响 Detached 会话继续写入。
    owner_core.close()
    assert rpcs["owner"].teardown_calls == [TEARDOWN_TIMEOUT]
    assert rpcs["detached"].teardown_calls == []

    detached_core.upsert_scalars([ScalarRecord(key="loss", step=1)])
    assert len(rpcs["detached"].calls) == 1

    detached_core.close()
    assert rpcs["detached"].teardown_calls == []


def test_core_go_rejects_negative_idle_timeout():
    with pytest.raises(ValueError):
        CoreGo("online", detach=True, idle_timeout=-1.0)


@pytest.mark.parametrize("reason", [ReasonCode.FAILED_PRECONDITION, ReasonCode.GRPC_UNAVAILABLE])
def test_detached_reuse_dying_core_spinup_retries_and_recovers(install_fake_launcher, monkeypatch, tmp_path, reason):
    """复用正在空闲退出的 Core 导致 spinup 失败时，应关闭旧通道并重试拉起。"""
    binary = tmp_path / "swanlab-core"
    attempts, _ = install_fake_launcher(
        core_module,
        owns_process=lambda index: index > 0,
        abort_increments_closed=True,
    )

    class FlakySpinupRpc(FakeRpc):
        def __init__(self, channel) -> None:
            super().__init__()
            self.channel = channel

        def spinup(self) -> None:
            if self.channel == "ch-0":
                raise CoreGoError(Stage.SPINUP, reason, "service is stopping or closed")

        def deliver_run_start(self, _req) -> DeliverRunStartResponse:
            return DeliverRunStartResponse(success=True, session_id="sid-fresh")

    monkeypatch.setattr(core_module, "CoreServiceRpc", FlakySpinupRpc)

    core = CoreGo("online", binary=binary, detach=True, idle_timeout=2.0)
    resp = core.deliver_run_start(DeliverRunStartRequest())

    assert resp.success is True
    assert len(attempts) == 2
    assert attempts[0].closed == 1
    core.close()


def test_owned_core_spinup_failure_does_not_retry(install_fake_launcher, monkeypatch, tmp_path):
    """自建进程 spinup 失败不应静默重试多次拉起。"""
    binary = tmp_path / "swanlab-core"
    attempts, starts = install_fake_launcher(core_module, owns_process=True, abort_increments_closed=True)

    class FailingRpc(FakeRpc):
        def __init__(self, _channel) -> None:
            super().__init__()

        def spinup(self) -> None:
            raise CoreGoError(Stage.SPINUP, ReasonCode.FAILED_PRECONDITION, "service is stopping or closed")

    monkeypatch.setattr(core_module, "CoreServiceRpc", FailingRpc)

    core = CoreGo("online", binary=binary, detach=True)
    resp = core.deliver_run_start(DeliverRunStartRequest())

    assert resp.success is False
    assert len(attempts) == 1
    assert len(starts) == 1


def test_detached_reuse_non_retryable_spinup_error_does_not_retry(install_fake_launcher, monkeypatch, tmp_path):
    """旧版本 Core 不支持 Spinup 时重试无意义，只报一次启动失败。"""
    binary = tmp_path / "swanlab-core"
    attempts, _ = install_fake_launcher(core_module, owns_process=False, abort_increments_closed=True)

    class UnimplementedSpinupRpc(FakeRpc):
        def __init__(self, _channel) -> None:
            super().__init__()

        def spinup(self) -> None:
            raise CoreGoError(Stage.SPINUP, ReasonCode.UNIMPLEMENTED, "spinup is not implemented")

    monkeypatch.setattr(core_module, "CoreServiceRpc", UnimplementedSpinupRpc)

    core = CoreGo("online", binary=binary, detach=True, idle_timeout=2.0)
    resp = core.deliver_run_start(DeliverRunStartRequest())

    assert resp.success is False
    assert len(attempts) == 1
    core.close()


@pytest.mark.parametrize(
    "stage, reason, grpc_name, expected_kind",
    [
        (Stage.RUN_START, ReasonCode.GRPC_UNAVAILABLE, "UNAVAILABLE", CORE_ERROR_KIND_TRANSIENT),
        (Stage.RUN_START, ReasonCode.GRPC_ERROR, "UNAUTHENTICATED", CORE_ERROR_KIND_AUTH),
        (Stage.RUN_START, ReasonCode.INVALID_ARGUMENT, "INVALID_ARGUMENT", CORE_ERROR_KIND_OTHER),
        (Stage.SPINUP, ReasonCode.GRPC_UNAVAILABLE, "UNAVAILABLE", CORE_ERROR_KIND_OTHER),
    ],
)
def test_deliver_run_start_populates_error_kind(
    install_fake_launcher, monkeypatch, tmp_path, stage, reason, grpc_name, expected_kind
):
    binary = tmp_path / "swanlab-core"
    install_fake_launcher(core_module, owns_process=True, abort_increments_closed=True)
    raised = CoreGoError(stage, reason, "start failed", grpc_name=grpc_name)

    class StartErrorRpc(FakeRpc):
        def __init__(self, _channel) -> None:
            super().__init__()

        def spinup(self) -> None:
            if stage is Stage.SPINUP:
                raise raised

        def deliver_run_start(self, _req) -> DeliverRunStartResponse:
            raise raised

    monkeypatch.setattr(core_module, "CoreServiceRpc", StartErrorRpc)

    core = CoreGo("online", binary=binary)
    resp = core.deliver_run_start(DeliverRunStartRequest())

    assert resp.success is False
    assert resp.error_kind == expected_kind
    assert core.last_error is raised


def test_deliver_run_start_lifecycle_failures_return_other_kind():
    closed_core = CoreGo("online")
    closed_core.close()
    closed_resp = closed_core.deliver_run_start(DeliverRunStartRequest())
    assert closed_resp.success is False
    assert closed_resp.error_kind == CORE_ERROR_KIND_OTHER

    running_core = _running_core(FakeRpc())
    started_resp = running_core.deliver_run_start(DeliverRunStartRequest())
    assert started_resp.success is False
    assert started_resp.error_kind == CORE_ERROR_KIND_OTHER


def test_confirm_run_finish_clears_stale_last_error():
    from swanlab.sdk.internal.run.finish_retry import retry_finish_call

    class ConfirmRpc(FakeRpc):
        def __init__(self, result: ConfirmRunFinishResponse, err: Optional[CoreGoError] = None) -> None:
            super().__init__(err)
            self.result = result
            self.confirm_calls = 0

        def confirm_run_finish(self, _session_id: str) -> ConfirmRunFinishResponse:
            self.confirm_calls += 1
            if self.error is not None:
                raise self.error
            return self.result

    # 1. 上一步遗留可重试错误，confirm 返回 success=False（未抛异常）时不应残留旧错误或触发重试
    stale = CoreGoError(Stage.RPC, ReasonCode.DEADLINE_EXCEEDED, "stats timeout")
    rpc_fail = ConfirmRpc(ConfirmRunFinishResponse(success=False, message="attach failed"))
    core = _running_core(rpc_fail)
    core._state = CoreState.FINISHED  # noqa: SLF001
    core._last_error = stale  # noqa: SLF001
    sleeps: List[float] = []
    resp = retry_finish_call(cast(CoreProtocol, core), core.confirm_run_finish, sleep=sleeps.append)
    assert resp.success is False
    assert resp.message == "attach failed"
    assert core.last_error is None
    assert rpc_fail.confirm_calls == 1
    assert sleeps == []

    # 2. confirm 成功也必须清零遗留错误
    rpc_ok = ConfirmRpc(ConfirmRunFinishResponse(success=True, message="OK"))
    core_ok = _running_core(rpc_ok)
    core_ok._state = CoreState.FINISHED  # noqa: SLF001
    core_ok._last_error = stale  # noqa: SLF001
    resp_ok = core_ok.confirm_run_finish()
    assert resp_ok.success is True
    assert core_ok.last_error is None

    # 3. confirm 自身抛出新的可重试错误时记录新错误
    fresh = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "transient confirm")
    rpc_err = ConfirmRpc(ConfirmRunFinishResponse(success=True), err=fresh)
    core_err = _running_core(rpc_err)
    core_err._state = CoreState.FINISHED  # noqa: SLF001
    core_err._last_error = stale  # noqa: SLF001
    resp_err = core_err.confirm_run_finish()
    assert resp_err.success is False
    assert core_err.last_error is fresh


class FakeFinishRpc:
    """只实现 finish / confirm / teardown，允许按需注入异常或响应。"""

    def __init__(self) -> None:
        self.deliver_error: Optional[BaseException] = None
        self.confirm_error: Optional[BaseException] = None
        self.deliver_calls = 0
        self.confirm_calls = 0
        self.teardown_calls: List[float] = []

    def deliver_run_finish(self, _request: DeliverRunFinishRequest, _session_id: str) -> DeliverRunFinishResponse:
        self.deliver_calls += 1
        if self.deliver_error is not None:
            raise self.deliver_error
        return DeliverRunFinishResponse(success=True)

    def confirm_run_finish(self, _session_id: str) -> ConfirmRunFinishResponse:
        self.confirm_calls += 1
        if self.confirm_error is not None:
            raise self.confirm_error
        return ConfirmRunFinishResponse(success=True)

    def teardown(self, timeout: float = 0) -> None:
        self.teardown_calls.append(timeout)


def _finished_core(rpc: FakeFinishRpc) -> CoreGo:
    core = _running_core(cast(FakeRpc, rpc))
    core._state = CoreState.FINISHED  # noqa: SLF001
    return core


def test_deliver_finish_rpc_error_is_retryable_and_keeps_running_state():
    error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "core unreachable")
    rpc = FakeFinishRpc()
    rpc.deliver_error = error
    core = _running_core(cast(FakeRpc, rpc))

    response = core.deliver_run_finish(DeliverRunFinishRequest())

    assert response.success is False
    assert core.last_error is error
    assert error.retryable is True
    assert core._state is CoreState.RUNNING  # noqa: SLF001


def test_deliver_finish_rejection_is_terminal():
    error = CoreGoError(Stage.RPC, ReasonCode.RUN_FINISH_REJECTED, "finish record commit failed")
    rpc = FakeFinishRpc()
    rpc.deliver_error = error
    core = _running_core(cast(FakeRpc, rpc))

    response = core.deliver_run_finish(DeliverRunFinishRequest())

    assert response.success is False
    assert core.last_error is error
    assert error.retryable is False


def test_confirm_rpc_error_is_retryable_and_keeps_session():
    error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "stop experiment failed")
    rpc = FakeFinishRpc()
    rpc.confirm_error = error
    core = _finished_core(rpc)

    response = core.confirm_run_finish()

    assert response.success is False
    assert core.last_error is error
    assert error.retryable is True
    assert core._session_id == "session-1"  # noqa: SLF001


def test_confirm_success_releases_session():
    rpc = FakeFinishRpc()
    core = _finished_core(rpc)

    response = core.confirm_run_finish()

    assert response.success is True
    assert core._session_id is None  # noqa: SLF001


def test_confirm_retry_after_retryable_failure_succeeds():
    rpc = FakeFinishRpc()
    rpc.confirm_error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "503")
    core = _finished_core(rpc)

    assert core.confirm_run_finish().success is False
    rpc.confirm_error = None
    assert core.confirm_run_finish().success is True
    assert rpc.confirm_calls == 2
    assert core._session_id is None  # noqa: SLF001


def test_close_after_confirm_failure_tears_down():
    rpc = FakeFinishRpc()
    rpc.confirm_error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "503")
    core = _finished_core(rpc)
    core.confirm_run_finish()

    core.close()

    assert core._state is CoreState.CLOSED  # noqa: SLF001
    assert rpc.teardown_calls == [TEARDOWN_TIMEOUT]


def test_close_records_core_exit_code(fake_launcher_cls):
    core = _finished_core(FakeFinishRpc())
    core._launcher = cast(CoreLauncher, fake_launcher_cls(exit_code=0))  # noqa: SLF001

    core.close()

    assert core.exit_code == 0


def test_core_go_forwards_log_dir_to_launcher(install_fake_launcher, monkeypatch, tmp_path):
    launchers, _ = install_fake_launcher(core_module)

    class StubRpc:
        def __init__(self, _channel) -> None:
            pass

        def spinup(self) -> None:
            pass

        def deliver_run_start(self, _request) -> DeliverRunStartResponse:
            return DeliverRunStartResponse(success=True, session_id="s-1")

        def teardown(self, timeout: float = 0.0) -> None:
            _ = timeout

    monkeypatch.setattr(core_module, "CoreServiceRpc", StubRpc)
    log_dir = tmp_path / "debug"
    core = CoreGo("local", binary=tmp_path / "swanlab-core", log_dir=log_dir)

    resp = core.deliver_run_start(DeliverRunStartRequest())

    assert resp.success is True
    assert launchers[0].log_dir == log_dir
    core.abort()
