"""
@author: Nexisato
@file: test_sync_service.py
@time: 2026/11/18
@description: CoreSyncRpc 的请求构造、超时语义与协议错误映射测试。
"""

from typing import Optional

import pytest
from google.protobuf.empty_pb2 import Empty

import swanlab.proto.swanlab.grpc.core.v1.sync_pb2_grpc as sync_pb2_grpc
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import GetOperationStatsResponse
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    ConfirmSyncFinishResponse,
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
    DeliverSyncStartResponse,
)
from swanlab.sdk.internal.core_go.config import (
    CONFIRM_TIMEOUT,
    RPC_TIMEOUT,
    RUN_START_TIMEOUT,
    STATS_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.rpc import CoreSyncRpc


@pytest.fixture
def sync_service(monkeypatch, fake_rpc_stub, fake_rpc_channel):
    monkeypatch.setattr(sync_pb2_grpc, "CoreSyncServiceStub", lambda raw: fake_rpc_stub)
    return CoreSyncRpc(fake_rpc_channel), fake_rpc_stub, fake_rpc_channel


def test_sync_service_module_does_not_load_grpc(run_python_isolated):
    code = (
        "import sys; import swanlab.sdk.internal.core_go.rpc.sync_service; "
        "print('grpc' in sys.modules, 'sync_pb2_grpc' in ' '.join(sys.modules))"
    )
    result = run_python_isolated(code)

    assert result.stdout.strip() == "False False", result.stderr


def test_deliver_sync_start_calls_stub_and_checks_success(sync_service):
    client, stub, channel = sync_service
    stub.responses["DeliverSyncStart"] = DeliverSyncStartResponse(success=True, message="success")
    req = DeliverSyncStartRequest(workspace="ws", project="proj", id="r1")

    resp = client.deliver_sync_start(req)

    assert resp.success is True
    name, sent_req, timeout = stub.calls[0]
    assert name == "DeliverSyncStart"
    assert sent_req is req
    assert timeout == RUN_START_TIMEOUT
    assert channel.stages[0] is Stage.RUN_START


def test_deliver_sync_start_rejected_raises_non_retryable_error(sync_service):
    client, stub, _ = sync_service
    stub.responses["DeliverSyncStart"] = DeliverSyncStartResponse(success=False, message="missing start record")

    with pytest.raises(CoreGoError) as caught:
        client.deliver_sync_start(DeliverSyncStartRequest())

    assert caught.value.stage is Stage.RUN_START
    assert caught.value.reason is ReasonCode.RUN_START_REJECTED
    assert caught.value.retryable is False
    assert "missing start record" in caught.value.message


def test_deliver_sync_flush_sends_empty_and_checks_success(sync_service):
    client, stub, channel = sync_service
    stub.responses["DeliverSyncFlush"] = DeliverSyncFlushResponse(success=True, path="/@alice/proj/slug")

    resp = client.deliver_sync_flush()

    assert resp.path == "/@alice/proj/slug"
    name, sent_req, timeout = stub.calls[0]
    assert name == "DeliverSyncFlush"
    assert isinstance(sent_req, Empty)
    assert timeout == RUN_START_TIMEOUT
    assert channel.stages[0] is Stage.RUN_START


def test_deliver_sync_flush_rejected_raises_error(sync_service):
    client, stub, _ = sync_service
    stub.responses["DeliverSyncFlush"] = DeliverSyncFlushResponse(success=False, message="403 forbidden")

    with pytest.raises(CoreGoError) as caught:
        client.deliver_sync_flush()

    assert caught.value.reason is ReasonCode.RUN_START_REJECTED
    assert caught.value.retryable is False
    assert "403 forbidden" in caught.value.message


def test_get_operation_stats_rejected_raises_protocol_violation(sync_service):
    client, stub, _ = sync_service
    stub.responses["GetOperationStats"] = GetOperationStatsResponse(success=False, message="sync session not started")

    with pytest.raises(CoreGoError) as caught:
        client.get_operation_stats()

    assert caught.value.reason is ReasonCode.PROTOCOL_VIOLATION


def test_confirm_sync_finish_rejected_raises_run_finish_rejected(sync_service):
    client, stub, _ = sync_service
    stub.responses["ConfirmSyncFinish"] = ConfirmSyncFinishResponse(success=False, message="stop failed")

    with pytest.raises(CoreGoError) as caught:
        client.confirm_sync_finish()

    assert caught.value.reason is ReasonCode.RUN_FINISH_REJECTED
    assert caught.value.retryable is False


def test_sync_rpc_timeouts_match_semantics(sync_service):
    client, stub, _ = sync_service
    stub.responses["DeliverSyncStart"] = DeliverSyncStartResponse(success=True)
    stub.responses["DeliverSyncFlush"] = DeliverSyncFlushResponse(success=True, path="/p")
    stub.responses["GetOperationStats"] = GetOperationStatsResponse(success=True)
    stub.responses["ConfirmSyncFinish"] = ConfirmSyncFinishResponse(success=True)

    client.deliver_sync_start(DeliverSyncStartRequest())
    client.deliver_sync_flush()
    client.get_operation_stats()
    client.confirm_sync_finish()

    timeouts = {name: timeout for name, _, timeout in stub.calls}
    assert timeouts["DeliverSyncStart"] == RUN_START_TIMEOUT
    assert timeouts["DeliverSyncFlush"] == RUN_START_TIMEOUT
    assert timeouts["GetOperationStats"] == STATS_TIMEOUT
    assert timeouts["ConfirmSyncFinish"] is CONFIRM_TIMEOUT is None


def test_unimplemented_is_not_rewritten(monkeypatch, fake_rpc_stub):
    class ErrorChannel:
        raw = object()

        def call(self, method, request, timeout: Optional[float] = RPC_TIMEOUT, *, stage=Stage.RPC):
            raise CoreGoError(
                Stage.RPC,
                ReasonCode.UNIMPLEMENTED,
                "CoreSyncService not implemented",
                grpc_name="UNIMPLEMENTED",
                grpc_code=12,
            )

    monkeypatch.setattr(sync_pb2_grpc, "CoreSyncServiceStub", lambda raw: fake_rpc_stub)
    client = CoreSyncRpc(ErrorChannel())

    with pytest.raises(CoreGoError) as caught:
        client.deliver_sync_flush()

    assert caught.value.reason is ReasonCode.UNIMPLEMENTED
    assert caught.value.grpc_name == "UNIMPLEMENTED"
    assert caught.value.grpc_code == 12
