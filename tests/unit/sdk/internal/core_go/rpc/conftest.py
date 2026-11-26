"""core_go RPC 单元测试共享 Stub / Channel fixture。"""

from typing import Optional

import pytest

from swanlab.sdk.internal.core_go.config import RPC_TIMEOUT, Stage


class FakeRpcStub:
    """按方法名记录调用并按需返回预置响应。"""

    def __init__(self) -> None:
        self.calls = []
        self.responses = {}

    def __getattr__(self, name: str):
        if name.startswith("_"):
            raise AttributeError(name)

        def method(request, timeout=None):
            self.calls.append((name, request, timeout))
            return self.responses.get(name)

        return method


class FakeRpcChannel:
    """记录调用 Stage 并直接透传给 FakeRpcStub 方法的假 Channel。"""

    def __init__(self) -> None:
        self.raw = object()
        self.stages = []

    def call(self, method, request, timeout: Optional[float] = RPC_TIMEOUT, *, stage=Stage.RPC):
        self.stages.append(stage)
        return method(request, timeout=timeout)


@pytest.fixture
def fake_rpc_stub() -> FakeRpcStub:
    return FakeRpcStub()


@pytest.fixture
def fake_rpc_channel() -> FakeRpcChannel:
    return FakeRpcChannel()
