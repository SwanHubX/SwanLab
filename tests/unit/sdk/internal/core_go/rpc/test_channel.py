"""
@author: Nexisato
@file: test_channel.py
@time: 2026/10/27
@description: Channel 的延迟导入、ready 超时与 gRPC 状态映射测试。
"""

import contextlib
import os
import socket
import sys
from concurrent.futures import ThreadPoolExecutor

import pytest

from swanlab.sdk.internal.core_go.config import GRPC_MAX_MESSAGE_BYTES, CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.rpc import Channel
from swanlab.sdk.internal.core_go.runtime import Endpoint


def test_rpc_package_does_not_load_grpc(run_python_isolated):
    """rpc 包顶层只加载 protobuf 消息，不加载 grpcio。"""
    code = (
        "import sys; import swanlab.sdk.internal.core_go.rpc; "
        "print('grpc' in sys.modules, 'core_pb2_grpc' in ' '.join(sys.modules))"
    )
    result = run_python_isolated(code)

    assert result.stdout.strip() == "False False", result.stderr


def test_channel_missing_grpcio_is_reported_as_local_error(monkeypatch):
    """grpcio 不在环境里时给出明确错误，而不是 ImportError 冒到上层。"""
    monkeypatch.setitem(sys.modules, "grpc", None)
    endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=1)

    with pytest.raises(CoreGoError) as caught:
        Channel.open(endpoint, timeout=0.1)

    assert caught.value.stage is Stage.CHANNEL
    # 本机缺少依赖重试无用，因此不能归到可重试的 GRPC_UNAVAILABLE
    assert caught.value.reason is ReasonCode.GRPC_LOCAL_ERROR
    assert caught.value.retryable is False
    assert "grpcio is required" in caught.value.message


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@contextlib.contextmanager
def _abort_server(code, details="boom"):
    grpc = pytest.importorskip("grpc")

    class Handler(grpc.GenericRpcHandler):
        def service(self, handler_call_details):
            def abort(request, context):
                context.abort(code, details)

            return grpc.unary_unary_rpc_method_handler(abort)

    server = grpc.server(ThreadPoolExecutor(max_workers=2))
    server.add_generic_rpc_handlers((Handler(),))
    port = server.add_insecure_port("127.0.0.1:0")
    server.start()
    try:
        yield port
    finally:
        server.stop(0)


def test_channel_ready_timeout():
    pytest.importorskip("grpc")
    endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=_free_port())

    with pytest.raises(CoreGoError) as caught:
        Channel.open(endpoint, timeout=0.2)

    assert caught.value.stage is Stage.CHANNEL
    assert caught.value.reason is ReasonCode.CHANNEL_TIMEOUT


def test_channel_open_disables_grpc_idle_timeout(monkeypatch):
    grpc = pytest.importorskip("grpc")
    captured = {}

    class Raw:
        def close(self):
            pass

    class Ready:
        def result(self, timeout=None):
            return None

    def fake_insecure_channel(target, options):
        captured["options"] = options
        return Raw()

    monkeypatch.setattr(grpc, "insecure_channel", fake_insecure_channel)
    monkeypatch.setattr(grpc, "channel_ready_future", lambda raw: Ready())
    endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=1)

    channel = Channel.open(endpoint, timeout=0.1)
    channel.close()

    assert ("grpc.client_idle_timeout_ms", 2**31 - 1) in captured["options"]
    assert ("grpc.enable_http_proxy", 0) in captured["options"]


def test_channel_open_ignores_http_proxy_env_on_loopback(monkeypatch):
    grpc = pytest.importorskip("grpc")
    monkeypatch.setenv("http_proxy", "http://127.0.0.1:1")
    monkeypatch.setenv("HTTP_PROXY", "http://127.0.0.1:1")
    monkeypatch.delenv("no_proxy", raising=False)
    monkeypatch.delenv("NO_PROXY", raising=False)

    with _abort_server(grpc.StatusCode.OK) as port:
        channel = Channel.open(Endpoint(kind="tcp", pid=os.getpid(), tcp_port=port), timeout=2.0)
        channel.close()


def test_channel_rejects_oversized_request_before_send():
    grpc = pytest.importorskip("grpc")
    channel = Channel(object(), grpc, "test")
    called = []

    class Oversized:
        def ByteSize(self):
            return GRPC_MAX_MESSAGE_BYTES + 1

    with pytest.raises(CoreGoError) as caught:
        channel.call(lambda request, timeout: called.append(request), Oversized(), timeout=0.1)

    assert called == []
    assert caught.value.reason is ReasonCode.REQUEST_TOO_LARGE
    assert caught.value.retryable is False
    assert caught.value.too_large is True


@pytest.mark.parametrize(
    "grpc_name, expected_reason, expected_code",
    [
        ("UNAVAILABLE", ReasonCode.GRPC_UNAVAILABLE, 14),
        ("UNIMPLEMENTED", ReasonCode.UNIMPLEMENTED, 12),
        ("NOT_FOUND", ReasonCode.SESSION_NOT_FOUND, 5),
    ],
)
def test_channel_maps_grpc_status(grpc_name, expected_reason, expected_code):
    grpc = pytest.importorskip("grpc")
    status = getattr(grpc.StatusCode, grpc_name)

    with _abort_server(status) as port:
        endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=port)
        channel = Channel.open(endpoint, timeout=2.0)
        method = channel.raw.unary_unary(
            "/smoke/abort",
            request_serializer=lambda value: value,
            response_deserializer=lambda value: value,
        )
        try:
            with pytest.raises(CoreGoError) as caught:
                channel.call(method, b"", timeout=1.0)
        finally:
            channel.close()

    error = caught.value
    assert error.stage is Stage.RPC
    assert error.reason is expected_reason
    assert error.grpc_name == grpc_name
    assert error.grpc_code == expected_code
    assert error.__cause__ is not None


@pytest.mark.parametrize(
    ("details", "expected_reason", "expected_retryable"),
    [
        ("Sent message larger than max (4194305 vs. 4194304)", ReasonCode.REQUEST_TOO_LARGE, False),
        ("server is busy", ReasonCode.RESOURCE_EXHAUSTED, True),
    ],
)
def test_resource_exhausted_message_too_large_is_terminal(details, expected_reason, expected_retryable):
    grpc = pytest.importorskip("grpc")

    with _abort_server(grpc.StatusCode.RESOURCE_EXHAUSTED, details) as port:
        endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=port)
        channel = Channel.open(endpoint, timeout=2.0)
        method = channel.raw.unary_unary(
            "/smoke/abort",
            request_serializer=lambda value: value,
            response_deserializer=lambda value: value,
        )
        try:
            with pytest.raises(CoreGoError) as caught:
                channel.call(method, b"", timeout=1.0)
        finally:
            channel.close()

    assert caught.value.reason is expected_reason
    assert caught.value.retryable is expected_retryable
    assert caught.value.too_large is (expected_reason is ReasonCode.REQUEST_TOO_LARGE)


def test_channel_close_is_idempotent_and_rejects_calls():
    grpc = pytest.importorskip("grpc")

    with _abort_server(grpc.StatusCode.UNAVAILABLE) as port:
        endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=port)
        channel = Channel.open(endpoint, timeout=2.0)
        channel.close()
        channel.close()

        with pytest.raises(CoreGoError) as caught:
            channel.call(lambda request, timeout: None, b"", timeout=1.0)

    assert caught.value.reason is ReasonCode.CLOSED


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [({}, "0"), ({"GRPC_ENABLE_FORK_SUPPORT": "1"}, "1")],
    ids=["default", "explicit"],
)
def test_fork_support_default_and_explicit_value(run_python_isolated, overrides, expected):
    """首次导入 rpc 包时设置默认值，且不覆盖用户显式配置。"""
    code = "import os; import swanlab.sdk.internal.core_go.rpc; print(os.environ['GRPC_ENABLE_FORK_SUPPORT'])"

    result = run_python_isolated(
        code,
        env_overrides=overrides,
        drop_env=("GRPC_ENABLE_FORK_SUPPORT",),
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == expected


@pytest.mark.skipif(sys.platform == "win32", reason="fork is POSIX only")
def test_fork_with_active_grpc_does_not_log_atfork_noise(run_python_isolated):
    """打开 fork 支持时，fork 会向 stderr 打印 absl 日志；默认值必须避免该噪音。"""
    code = "\n".join(
        [
            "import os, threading, time",
            "import swanlab.sdk.internal.core_go.rpc",
            "import grpc",
            "channel = grpc.insecure_channel('10.255.255.1:1')",
            "def probe():",
            "    try:",
            "        grpc.channel_ready_future(channel).result(timeout=3)",
            "    except Exception:",
            "        pass",
            "threading.Thread(target=probe, daemon=True).start()",
            "time.sleep(0.3)",
            "pid = os.fork()",
            "if pid == 0:",
            "    os._exit(0)",
            "os.waitpid(pid, 0)",
        ]
    )

    result = run_python_isolated(code, drop_env=("GRPC_ENABLE_FORK_SUPPORT",), check=False)

    assert result.returncode == 0, result.stderr
    assert "fork_posix.cc" not in result.stderr
