"""
@author: Nexisato
@file: test_channel.py
@time: 2026/10/27
@description: Channel 的延迟导入、ready 超时与 gRPC 状态映射测试。
"""

import contextlib
import os
import socket
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.rpc import Channel
from swanlab.sdk.internal.core_go.runtime import Endpoint

REPO_ROOT = Path(__file__).resolve().parents[5]


def test_rpc_package_does_not_load_grpc():
    """rpc 包顶层只加载 protobuf 消息，不加载 grpcio。"""
    code = (
        "import sys; import swanlab.sdk.internal.core_go.rpc; "
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


def test_channel_missing_grpcio_is_reported_as_unavailable(monkeypatch):
    """grpcio 不在环境里时给出明确错误，而不是 ImportError 冒到上层。"""
    monkeypatch.setitem(sys.modules, "grpc", None)
    endpoint = Endpoint(kind="tcp", pid=os.getpid(), tcp_port=1)

    with pytest.raises(CoreGoError) as caught:
        Channel.open(endpoint, timeout=0.1)

    assert caught.value.stage is Stage.CHANNEL
    assert caught.value.reason is ReasonCode.GRPC_UNAVAILABLE
    assert "grpcio is required" in caught.value.message


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@contextlib.contextmanager
def _abort_server(code):
    grpc = pytest.importorskip("grpc")

    class Handler(grpc.GenericRpcHandler):
        def service(self, handler_call_details):
            def abort(request, context):
                context.abort(code, "boom")

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


def _run_python(code: str, *, env_overrides: dict) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k != "GRPC_ENABLE_FORK_SUPPORT"}
    env["PYTHONPATH"] = str(REPO_ROOT)
    env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
    )


def test_fork_support_default_is_set_before_grpc_import():
    """首次导入 rpc 包时关闭 grpc fork 支持，避免 fork 子进程 stderr 噪音。"""
    code = "import os; import swanlab.sdk.internal.core_go.rpc; print(os.environ['GRPC_ENABLE_FORK_SUPPORT'])"

    result = _run_python(code, env_overrides={})

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "0"


def test_fork_support_explicit_value_is_preserved():
    code = "import os; import swanlab.sdk.internal.core_go.rpc; print(os.environ['GRPC_ENABLE_FORK_SUPPORT'])"

    result = _run_python(code, env_overrides={"GRPC_ENABLE_FORK_SUPPORT": "1"})

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.skipif(sys.platform == "win32", reason="fork is POSIX only")
def test_fork_with_active_grpc_does_not_log_atfork_noise():
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

    result = _run_python(code, env_overrides={})

    assert result.returncode == 0, result.stderr
    assert "fork_posix.cc" not in result.stderr
