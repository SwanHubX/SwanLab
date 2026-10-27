"""
@author: Nexisato
@file: test_core_service.py
@time: 2026/10/27
@description: CoreServiceRpc 的请求构造、session_id 与协议错误测试。
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

import swanlab.proto.swanlab.grpc.core.v1.core_pb2_grpc as core_pb2_grpc
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    ConfirmRunFinishResponse,
    DeliverRunFinishRequest,
    DeliverRunFinishResponse,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
    GetOperationStatsResponse,
    SpinupServiceRequest,
    UpsertColumnsRequest,
    UpsertLogsRequest,
    UpsertMediaRequest,
    UpsertSavesRequest,
    UpsertScalarsRequest,
)
from swanlab.sdk.internal.core_go.config import RPC_TIMEOUT
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.rpc import CoreServiceRpc

REPO_ROOT = Path(__file__).resolve().parents[5]


class FakeStub:
    """按方法名记录调用并按需返回预置响应。"""

    def __init__(self):
        self.calls = []
        self.responses = {}

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)

        def method(request, timeout=None):
            self.calls.append((name, request, timeout))
            return self.responses.get(name)

        return method


class FakeChannel:
    def __init__(self):
        self.raw = object()

    def call(self, method, request, timeout=RPC_TIMEOUT, *, stage=Stage.RPC):
        return method(request, timeout=timeout)


@pytest.fixture
def service(monkeypatch):
    stub = FakeStub()
    monkeypatch.setattr(core_pb2_grpc, "CoreServiceStub", lambda raw: stub)
    return CoreServiceRpc(FakeChannel()), stub


def test_core_service_module_does_not_load_grpc():
    code = (
        "import sys; import swanlab.sdk.internal.core_go.rpc.core_service; "
        "print('grpc' in sys.modules, 'core_pb2_grpc' in ' '.join(sys.modules))"
    )
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )

    assert result.stdout.strip() == "False False", result.stderr


def test_spinup_calls_stub_once(service):
    client, stub = service

    client.spinup()

    assert [call[0] for call in stub.calls] == ["SpinupService"]
    assert isinstance(stub.calls[0][1], SpinupServiceRequest)


def test_teardown_uses_given_timeout(service):
    client, stub = service

    client.teardown(timeout=1.5)

    name, request, timeout = stub.calls[0]
    assert name == "TeardownService"
    assert timeout == 1.5
    assert request is not None


def test_deliver_run_start_returns_session(service):
    client, stub = service
    stub.responses["DeliverRunStart"] = DeliverRunStartResponse(success=True, session_id="s-1")

    response = client.deliver_run_start(DeliverRunStartRequest())

    assert response.session_id == "s-1"
    assert [call[0] for call in stub.calls] == ["DeliverRunStart"]


def test_deliver_run_start_rejected(service):
    client, stub = service
    stub.responses["DeliverRunStart"] = DeliverRunStartResponse(success=False, message="bad run")

    with pytest.raises(CoreGoError) as caught:
        client.deliver_run_start(DeliverRunStartRequest())

    assert caught.value.stage is Stage.RUN_START
    assert caught.value.reason is ReasonCode.RUN_START_REJECTED
    assert caught.value.message == "bad run"


def test_deliver_run_start_requires_session(service):
    client, stub = service
    stub.responses["DeliverRunStart"] = DeliverRunStartResponse(success=True, session_id="")

    with pytest.raises(CoreGoError) as caught:
        client.deliver_run_start(DeliverRunStartRequest())

    assert caught.value.reason is ReasonCode.PROTOCOL_VIOLATION


@pytest.mark.parametrize(
    "method_name, stub_name, request_type, field",
    [
        ("upsert_columns", "UpsertColumns", UpsertColumnsRequest, "columns"),
        ("upsert_scalars", "UpsertScalars", UpsertScalarsRequest, "data"),
        ("upsert_media", "UpsertMedia", UpsertMediaRequest, "data"),
        ("upsert_logs", "UpsertLogs", UpsertLogsRequest, "logs"),
        ("upsert_saves", "UpsertSaves", UpsertSavesRequest, "saves"),
    ],
)
def test_upsert_attaches_session_id(service, method_name, stub_name, request_type, field):
    client, stub = service

    getattr(client, method_name)([], "s-1")

    name, request, _ = stub.calls[0]
    assert name == stub_name
    assert isinstance(request, request_type)
    assert request.session_id == "s-1"
    assert list(getattr(request, field)) == []


def test_deliver_run_finish_sets_session_id(service):
    client, stub = service
    stub.responses["DeliverRunFinish"] = DeliverRunFinishResponse(success=True)
    request = DeliverRunFinishRequest()

    client.deliver_run_finish(request, "s-1")

    assert request.session_id == "s-1"


def test_deliver_run_finish_rejected(service):
    client, stub = service
    stub.responses["DeliverRunFinish"] = DeliverRunFinishResponse(success=False, message="still running")

    with pytest.raises(CoreGoError) as caught:
        client.deliver_run_finish(DeliverRunFinishRequest(), "s-1")

    assert caught.value.reason is ReasonCode.RUN_FINISH_REJECTED


def test_get_operation_stats_rejected(service):
    client, stub = service
    stub.responses["GetOperationStats"] = GetOperationStatsResponse(success=False, message="no such session")

    with pytest.raises(CoreGoError) as caught:
        client.get_operation_stats("s-1")

    assert caught.value.reason is ReasonCode.PROTOCOL_VIOLATION


def test_confirm_run_finish_rejected(service):
    client, stub = service
    stub.responses["ConfirmRunFinish"] = ConfirmRunFinishResponse(success=False, message="not finished")

    with pytest.raises(CoreGoError) as caught:
        client.confirm_run_finish("s-1")

    assert caught.value.reason is ReasonCode.RUN_FINISH_REJECTED


def test_unimplemented_is_not_rewritten(monkeypatch):
    class ErrorChannel:
        raw = object()

        def call(self, method, request, timeout=RPC_TIMEOUT, *, stage=Stage.RPC):
            raise CoreGoError(
                Stage.RPC,
                ReasonCode.UNIMPLEMENTED,
                "UpsertScalars is not implemented yet",
                grpc_name="UNIMPLEMENTED",
                grpc_code=12,
            )

    monkeypatch.setattr(core_pb2_grpc, "CoreServiceStub", lambda raw: FakeStub())
    client = CoreServiceRpc(ErrorChannel())

    with pytest.raises(CoreGoError) as caught:
        client.upsert_scalars([], "s-1")

    assert caught.value.reason is ReasonCode.UNIMPLEMENTED
    assert caught.value.grpc_name == "UNIMPLEMENTED"
    assert caught.value.grpc_code == 12
