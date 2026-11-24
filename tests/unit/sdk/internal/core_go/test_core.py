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
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.core_go import core as core_module
from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.core import CoreGo, _CoreState, _DetachedCoreGo
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc
from swanlab.sdk.internal.pkg import console as pkg_console
from swanlab.sdk.protocol import CoreProtocol, CoreProtocolError, NullCoreProtocol, is_terminal_stats_error


class FakeRpc:
    """只实现 upsert_scalars，按需抛出预置异常。"""

    def __init__(self, error: Optional[BaseException] = None) -> None:
        self.error = error
        self.calls: List[tuple] = []
        self.teardown_calls: List[float] = []

    def upsert_scalars(self, scalars: List[ScalarRecord], session_id: str) -> None:
        self.calls.append((scalars, session_id))
        if self.error is not None:
            raise self.error

    def teardown(self, timeout: float = 0) -> None:
        self.teardown_calls.append(timeout)


def _running_core(rpc: FakeRpc) -> CoreGo:
    core = CoreGo("online")
    core._rpc = cast(CoreServiceRpc, rpc)  # noqa: SLF001 - 测试直接注入 fake RPC
    core._session_id = "session-1"  # noqa: SLF001
    core._state = _CoreState.RUNNING  # noqa: SLF001
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


def test_unknown_exception_is_wrapped_as_retryable():
    core = _running_core(FakeRpc(RuntimeError("socket closed")))

    with pytest.raises(CoreGoError) as caught:
        core.upsert_scalars([ScalarRecord(key="loss", step=1)])

    assert caught.value.reason is ReasonCode.GRPC_ERROR
    assert caught.value.retryable is True
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
    core._state = _CoreState.CLOSED  # noqa: SLF001

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

    assert core._state is _CoreState.CLOSED  # noqa: SLF001


def test_close_after_finish_tears_down_and_releases_rpc():
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._state = _CoreState.FINISHED  # noqa: SLF001

    core.close()

    assert core._state is _CoreState.CLOSED  # noqa: SLF001
    assert core._rpc is None  # noqa: SLF001
    assert rpc.teardown_calls == [core_module.TEARDOWN_TIMEOUT]


def test_fork_returns_detached_placeholder(tmp_path):
    binary = tmp_path / "swanlab-core"
    core = CoreGo("online", binary=binary)

    detached = core.fork()

    assert isinstance(detached, _DetachedCoreGo)
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
    monkeypatch.setattr(CoreGo, "_start_enabled", lambda self, req: calls.append(self._mode) or sentinel)

    for mode in ("online", "local", "offline"):
        assert CoreGo(mode).deliver_run_start(DeliverRunStartRequest()) is sentinel
    assert calls == ["online", "local", "offline"]


def test_finish_modes_share_enabled_path(monkeypatch):
    sentinel = DeliverRunFinishResponse(success=True)
    monkeypatch.setattr(CoreGo, "_finish_enabled", lambda self, req: sentinel)

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


class _SpyLauncher:
    def __init__(self, *, owns_process: bool = True) -> None:
        self.owns_process = owns_process
        self.closed = 0
        self.aborted = 0
        self.exit_code = None

    def close(self) -> None:
        self.closed += 1

    def abort(self) -> None:
        self.aborted += 1


@pytest.mark.parametrize("action, initial_state", [("close", _CoreState.FINISHED), ("abort", _CoreState.RUNNING)])
def test_detached_mode_close_and_abort_skip_teardown_and_are_idempotent(action, initial_state):
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._detach = True  # noqa: SLF001
    launcher = _SpyLauncher()
    core._launcher = cast(core_module.CoreLauncher, launcher)  # noqa: SLF001
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
    assert core._state is _CoreState.CLOSED  # noqa: SLF001


def test_owner_mode_abort_sends_teardown():
    rpc = FakeRpc()
    core = _running_core(rpc)

    core.abort(RuntimeError("boom"))

    assert rpc.teardown_calls == [core_module.TEARDOWN_TIMEOUT]
    assert core._state is _CoreState.CLOSED  # noqa: SLF001


def test_unowned_launcher_skips_teardown_on_close_even_when_detach_false():
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._detach = False  # noqa: SLF001
    launcher = _SpyLauncher(owns_process=False)
    core._launcher = cast(core_module.CoreLauncher, launcher)  # noqa: SLF001

    core.close()

    assert rpc.teardown_calls == []
    assert launcher.closed == 1
    assert core._state is _CoreState.CLOSED  # noqa: SLF001


@pytest.mark.parametrize("first_mode", ["detached", "owner"])
def test_mixed_owner_and_detached_close_does_not_affect_peer(monkeypatch, tmp_path, first_mode):
    binary = tmp_path / "swanlab-core"
    rpcs = {}
    launcher_kwargs = {}

    class StubLauncher:
        def __init__(self, bin_path, *, detach: bool = False, idle_timeout: float = 0.0) -> None:
            launcher_kwargs[detach] = {"binary": bin_path, "idle_timeout": idle_timeout}
            self.detach = detach
            self.owns_process = True
            self.channel = "detached" if detach else "owner"
            self.exit_code = None if detach else 0
            self.closed = 0

        def start(self) -> None:
            return None

        def close(self) -> None:
            self.closed += 1

        def abort(self) -> None:
            self.closed += 1

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

    monkeypatch.setattr(core_module, "CoreLauncher", StubLauncher)
    monkeypatch.setattr(core_module, "CoreServiceRpc", TaggedRpc)

    detached_core = CoreGo("online", binary=binary, detach=True, idle_timeout=2.0)
    owner_core = CoreGo("online", binary=binary, detach=False)

    if first_mode == "detached":
        assert detached_core.deliver_run_start(DeliverRunStartRequest()).success is True
        assert owner_core.deliver_run_start(DeliverRunStartRequest()).success is True
    else:
        assert owner_core.deliver_run_start(DeliverRunStartRequest()).success is True
        assert detached_core.deliver_run_start(DeliverRunStartRequest()).success is True

    assert launcher_kwargs == {
        True: {"binary": binary, "idle_timeout": 2.0},
        False: {"binary": binary, "idle_timeout": core_module.DETACH_IDLE_TIMEOUT},
    }

    # Owner 退出只 teardown 自己的 Core，不得影响 Detached 会话继续写入。
    owner_core.close()
    assert rpcs["owner"].teardown_calls == [core_module.TEARDOWN_TIMEOUT]
    assert rpcs["detached"].teardown_calls == []

    detached_core.upsert_scalars([ScalarRecord(key="loss", step=1)])
    assert len(rpcs["detached"].calls) == 1

    detached_core.close()
    assert rpcs["detached"].teardown_calls == []


def test_core_go_rejects_negative_idle_timeout():
    with pytest.raises(ValueError):
        CoreGo("online", detach=True, idle_timeout=-1.0)


def _install_stub_launcher(monkeypatch, *, owns_process):
    """安装共享的假 Launcher，返回 (attempts, starts)。

    ``owns_process`` 可为布尔值或 ``index -> bool``，用于覆盖复用与自建进程两条路径。
    """
    attempts = []
    starts = []

    class StubLauncher:
        def __init__(self, bin_path, *, detach: bool = False, idle_timeout: float = 0.0) -> None:
            self.index = len(attempts)
            self.owns_process = owns_process(self.index) if callable(owns_process) else owns_process
            self.channel = f"ch-{self.index}"
            self.exit_code = None
            self.closed = 0
            attempts.append(self)

        def start(self) -> None:
            starts.append(self.index)

        def close(self) -> None:
            self.closed += 1

        def abort(self) -> None:
            self.closed += 1

    monkeypatch.setattr(core_module, "CoreLauncher", StubLauncher)
    return attempts, starts


@pytest.mark.parametrize("reason", [ReasonCode.FAILED_PRECONDITION, ReasonCode.GRPC_UNAVAILABLE])
def test_detached_reuse_dying_core_spinup_retries_and_recovers(monkeypatch, tmp_path, reason):
    """复用正在空闲退出的 Core 导致 spinup 失败时，应关闭旧通道并重试拉起。"""
    binary = tmp_path / "swanlab-core"
    attempts, _ = _install_stub_launcher(monkeypatch, owns_process=lambda index: index > 0)

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


def test_owned_core_spinup_failure_does_not_retry(monkeypatch, tmp_path):
    """自建进程 spinup 失败不应静默重试多次拉起。"""
    binary = tmp_path / "swanlab-core"
    attempts, starts = _install_stub_launcher(monkeypatch, owns_process=True)

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


def test_detached_reuse_non_retryable_spinup_error_does_not_retry(monkeypatch, tmp_path):
    """旧版本 Core 不支持 Spinup 时重试无意义，只报一次启动失败。"""
    binary = tmp_path / "swanlab-core"
    attempts, _ = _install_stub_launcher(monkeypatch, owns_process=False)

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
        (Stage.RUN_START, ReasonCode.DEADLINE_EXCEEDED, "DEADLINE_EXCEEDED", CORE_ERROR_KIND_TRANSIENT),
        (Stage.RUN_START, ReasonCode.GRPC_ERROR, "UNAUTHENTICATED", CORE_ERROR_KIND_AUTH),
        (Stage.RUN_START, ReasonCode.INVALID_ARGUMENT, "INVALID_ARGUMENT", CORE_ERROR_KIND_OTHER),
        (Stage.RUN_START, ReasonCode.GRPC_ERROR, "INTERNAL", CORE_ERROR_KIND_OTHER),
        (Stage.RUN_START, ReasonCode.SESSION_NOT_FOUND, "NOT_FOUND", CORE_ERROR_KIND_OTHER),
        (Stage.RUN_START, ReasonCode.RUN_START_REJECTED, None, CORE_ERROR_KIND_OTHER),
        (Stage.SPINUP, ReasonCode.GRPC_UNAVAILABLE, "UNAVAILABLE", CORE_ERROR_KIND_OTHER),
        (Stage.CHANNEL, ReasonCode.DEADLINE_EXCEEDED, "DEADLINE_EXCEEDED", CORE_ERROR_KIND_OTHER),
    ],
)
def test_deliver_run_start_populates_error_kind(monkeypatch, tmp_path, stage, reason, grpc_name, expected_kind):
    binary = tmp_path / "swanlab-core"
    _install_stub_launcher(monkeypatch, owns_process=True)
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
    core._state = _CoreState.FINISHED  # noqa: SLF001
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
    core_ok._state = _CoreState.FINISHED  # noqa: SLF001
    core_ok._last_error = stale  # noqa: SLF001
    resp_ok = core_ok.confirm_run_finish()
    assert resp_ok.success is True
    assert core_ok.last_error is None

    # 3. confirm 自身抛出新的可重试错误时记录新错误
    fresh = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "transient confirm")
    rpc_err = ConfirmRpc(ConfirmRunFinishResponse(success=True), err=fresh)
    core_err = _running_core(rpc_err)
    core_err._state = _CoreState.FINISHED  # noqa: SLF001
    core_err._last_error = stale  # noqa: SLF001
    resp_err = core_err.confirm_run_finish()
    assert resp_err.success is False
    assert core_err.last_error is fresh
