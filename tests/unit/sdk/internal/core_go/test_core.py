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
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
)
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.core_go import core as core_module
from swanlab.sdk.internal.core_go.core import CoreGo, _DetachedCoreGo
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.rpc.core_service import CoreServiceRpc
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
    core._state = core_module._RUNNING  # noqa: SLF001
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
    core._state = core_module._CLOSED  # noqa: SLF001

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

    assert core._state == core_module._CLOSED  # noqa: SLF001


def test_close_after_finish_tears_down_and_releases_rpc():
    rpc = FakeRpc()
    core = _running_core(rpc)
    core._state = core_module._FINISHED  # noqa: SLF001

    core.close()

    assert core._state == core_module._CLOSED  # noqa: SLF001
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


def test_disabled_start_is_noop_without_binary():
    core = CoreGo("disabled")

    response = core.deliver_run_start(DeliverRunStartRequest())

    assert response.success is True
