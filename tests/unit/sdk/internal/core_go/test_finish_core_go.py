"""finish / confirm 失败分类与可重试语义。

关键契约：
- deliver 的确定性拒绝（如终态提交失败）不重试；
- confirm 的释放失败可重试，且 CoreGo 必须保留 session_id；
- 只有 confirm 成功才释放本地会话引用。
"""

from typing import List, Optional, cast

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
)
from swanlab.sdk.internal.core_go import core as core_module
from swanlab.sdk.internal.core_go.core import CoreGo
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.launcher import CoreLauncher
from swanlab.sdk.internal.core_go.rpc.core_service import CoreServiceRpc


class FakeFinishRpc:
    """只实现 finish / confirm / teardown，允许按需注入异常或响应。"""

    def __init__(self) -> None:
        self.deliver_error: Optional[BaseException] = None
        self.deliver_response = DeliverRunFinishResponse(success=True)
        self.confirm_error: Optional[BaseException] = None
        self.confirm_response = ConfirmRunFinishResponse(success=True)
        self.deliver_calls = 0
        self.confirm_calls = 0
        self.teardown_calls: List[float] = []

    def deliver_run_finish(self, _request: DeliverRunFinishRequest, _session_id: str) -> DeliverRunFinishResponse:
        self.deliver_calls += 1
        if self.deliver_error is not None:
            raise self.deliver_error
        return self.deliver_response

    def confirm_run_finish(self, _session_id: str) -> ConfirmRunFinishResponse:
        self.confirm_calls += 1
        if self.confirm_error is not None:
            raise self.confirm_error
        return self.confirm_response

    def teardown(self, timeout: float = 0) -> None:
        self.teardown_calls.append(timeout)


def _running_core(rpc: FakeFinishRpc) -> CoreGo:
    core = CoreGo("online")
    core._rpc = cast(CoreServiceRpc, rpc)  # noqa: SLF001 - 测试直接注入 fake RPC
    core._session_id = "session-1"  # noqa: SLF001
    core._state = core_module._RUNNING  # noqa: SLF001
    return core


def _finished_core(rpc: FakeFinishRpc) -> CoreGo:
    core = _running_core(rpc)
    core._state = core_module._FINISHED  # noqa: SLF001
    return core


def test_deliver_finish_rpc_error_is_retryable_and_keeps_running_state():
    error = CoreGoError(Stage.RPC, ReasonCode.GRPC_UNAVAILABLE, "core unreachable")
    rpc = FakeFinishRpc()
    rpc.deliver_error = error
    core = _running_core(rpc)

    response = core.deliver_run_finish(DeliverRunFinishRequest())

    assert response.success is False
    assert core.last_error is error
    assert error.retryable is True
    assert core._state == core_module._RUNNING  # noqa: SLF001


def test_deliver_finish_rejection_is_terminal():
    error = CoreGoError(Stage.RPC, ReasonCode.RUN_FINISH_REJECTED, "finish record commit failed")
    rpc = FakeFinishRpc()
    rpc.deliver_error = error
    core = _running_core(rpc)

    response = core.deliver_run_finish(DeliverRunFinishRequest())

    assert response.success is False
    assert core.last_error is error
    assert error.retryable is False


def test_deliver_finish_rejected_response_is_terminal():
    """RPC 层不抛错时，CoreGo 仍需把语义拒绝标记为不可重试。"""
    rpc = FakeFinishRpc()
    rpc.deliver_response = DeliverRunFinishResponse(success=False, message="store failed")
    core = _running_core(rpc)

    response = core.deliver_run_finish(DeliverRunFinishRequest())

    assert response.success is False
    assert core.last_error is not None
    assert core.last_error.reason is ReasonCode.RUN_FINISH_REJECTED
    assert core.last_error.retryable is False


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


def test_confirm_rejected_response_is_retryable_and_keeps_session():
    rpc = FakeFinishRpc()
    rpc.confirm_response = ConfirmRunFinishResponse(success=False, message="saved locally")
    core = _finished_core(rpc)

    response = core.confirm_run_finish()

    assert response.success is False
    assert core.last_error is not None
    assert core.last_error.reason is ReasonCode.RUN_FINISH_REJECTED
    assert core.last_error.retryable is True
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

    assert core._state == core_module._CLOSED  # noqa: SLF001
    assert rpc.teardown_calls == [core_module.TEARDOWN_TIMEOUT]


def test_close_records_core_exit_code():
    class FakeLauncher:
        exit_code = 0

        def close(self) -> None:
            pass

    core = _finished_core(FakeFinishRpc())
    core._launcher = cast(CoreLauncher, FakeLauncher())  # noqa: SLF001

    core.close()

    assert core.exit_code == 0
