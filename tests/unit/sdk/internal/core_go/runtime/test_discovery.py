"""
@author: Nexisato
@file: test_discovery.py
@time: 2026/11/16
@description: 已有 port-file 发现、端点存活探测与边界校验测试。
"""

import socket
import sys

import pytest

from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.runtime import (
    Endpoint,
    discover_endpoint,
    probe_endpoint,
    validate_discovered_endpoint,
)


def test_discover_endpoint_returns_none_when_missing_or_corrupt(short_tmp_path):
    port_file = short_tmp_path / "core.port"
    assert discover_endpoint(port_file, short_tmp_path) is None

    port_file.write_bytes(b"corrupt-content")
    assert discover_endpoint(port_file, short_tmp_path) is None


def test_discover_endpoint_finds_live_tcp_and_ignores_stale_tcp(short_tmp_path):
    port_file = short_tmp_path / "core.port"
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind(("127.0.0.1", 0))
        listener.listen(1)
        port = listener.getsockname()[1]

        port_file.write_bytes(f"protocol=1\npid=4321\nsock={port}\nEOF\n".encode("utf-8"))
        found = discover_endpoint(port_file, short_tmp_path)
        assert found == Endpoint(kind="tcp", pid=4321, tcp_port=port)

    # 监听关闭后，陈旧 port-file 不再被当成存活实例。
    assert discover_endpoint(port_file, short_tmp_path) is None


@pytest.mark.skipif(sys.platform == "win32" or not hasattr(socket, "AF_UNIX"), reason="UDS is POSIX-only")
def test_discover_endpoint_finds_live_unix_and_rejects_outside_dir(short_tmp_path, tmp_path):
    port_file = short_tmp_path / "core.port"
    sock_path = short_tmp_path / "core.sock"

    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as listener:
        listener.bind(str(sock_path))
        listener.listen(1)

        port_file.write_bytes(f"protocol=1\npid=8888\nunix={sock_path}\nEOF\n".encode("utf-8"))
        found = discover_endpoint(port_file, short_tmp_path)
        assert found == Endpoint(kind="unix", pid=8888, unix_path=str(sock_path))

    assert probe_endpoint(Endpoint(kind="unix", pid=8888, unix_path=str(sock_path))) is False

    outside_sock = tmp_path / "outside.sock"
    with pytest.raises(CoreGoError) as caught:
        validate_discovered_endpoint(
            Endpoint(kind="unix", pid=8888, unix_path=str(outside_sock)),
            short_tmp_path,
        )
    assert caught.value.stage is Stage.PORT_FILE
    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH

    port_file.write_bytes(f"protocol=1\npid=8888\nunix={outside_sock}\nEOF\n".encode("utf-8"))
    assert discover_endpoint(port_file, short_tmp_path) is None
