"""
@author: Nexisato
@file: test_ownership.py
@time: 2026/10/27
@description: Endpoint 归属校验测试。
"""

import sys

import pytest

from swanlab.sdk.internal.core_go.config import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.core_go.runtime import Endpoint, validate_ownership


def test_accepts_matching_unix_endpoint(tmp_path):
    socket_path = tmp_path / "core.sock"
    socket_path.touch()

    validate_ownership(Endpoint(kind="unix", pid=123, unix_path=str(socket_path)), 123, tmp_path)


def test_accepts_matching_tcp_endpoint(tmp_path):
    validate_ownership(Endpoint(kind="tcp", pid=123, tcp_port=50051), 123, tmp_path)


def test_rejects_pid_mismatch(tmp_path):
    socket_path = tmp_path / "core.sock"
    socket_path.touch()

    with pytest.raises(CoreGoError) as caught:
        validate_ownership(Endpoint(kind="unix", pid=124, unix_path=str(socket_path)), 123, tmp_path)

    assert caught.value.stage is Stage.PORT_FILE
    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH


def test_rejects_unix_socket_outside_runtime_dir(tmp_path):
    other_dir = tmp_path / "other"
    other_dir.mkdir()
    socket_path = other_dir / "core.sock"
    socket_path.touch()

    with pytest.raises(CoreGoError) as caught:
        validate_ownership(Endpoint(kind="unix", pid=123, unix_path=str(socket_path)), 123, tmp_path / "runtime")

    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX symlink check")
def test_rejects_unix_socket_symlink(tmp_path):
    runtime_dir = tmp_path / "runtime"
    runtime_dir.mkdir()
    outside = tmp_path / "outside.sock"
    outside.touch()
    (runtime_dir / "core.sock").symlink_to(outside)

    with pytest.raises(CoreGoError) as caught:
        validate_ownership(
            Endpoint(kind="unix", pid=123, unix_path=str(runtime_dir / "core.sock")),
            123,
            runtime_dir,
        )

    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH


def test_rejects_invalid_tcp_port(tmp_path):
    with pytest.raises(CoreGoError) as caught:
        validate_ownership(Endpoint(kind="tcp", pid=123, tcp_port=0), 123, tmp_path)

    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH


def test_rejects_unknown_endpoint_kind(tmp_path):
    with pytest.raises(CoreGoError) as caught:
        validate_ownership(Endpoint(kind="pipe", pid=123), 123, tmp_path)

    assert caught.value.reason is ReasonCode.OWNERSHIP_MISMATCH
