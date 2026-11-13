"""
@author: Nexisato
@file: test_port_file.py
@time: 2026/10/27
@description: port-file v1 解析测试，用例对齐 Go portinfo_test.go。
"""

import pytest

from swanlab.sdk.internal.core_go.config import MAX_UNIX_PATH_BYTES, PORT_FILE_MAX_BYTES
from swanlab.sdk.internal.core_go.discovery import Endpoint, parse_port_file, read_port_file
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage

VALID_UNIX = "protocol=1\npid=4242\nunix=/tmp/swanlab/core.sock\nEOF\n"


def test_parse_accepts_valid_variants():
    """TestCase: TestParseAcceptsValidVariants"""
    # EOF 行末尾可无换行
    endpoint = parse_port_file(VALID_UNIX.rstrip("\n").encode())
    assert endpoint == Endpoint(kind="unix", pid=4242, unix_path="/tmp/swanlab/core.sock")

    # 键值行顺序不影响解析
    reordered = b"unix=/tmp/swanlab/core.sock\npid=4242\nprotocol=1\nEOF\n"
    assert parse_port_file(reordered).unix_path == "/tmp/swanlab/core.sock"

    # sock 端点（Windows 或 UDS 回退形态）
    sock_content = b"protocol=1\npid=4242\nsock=12345\nEOF\n"
    assert parse_port_file(sock_content) == Endpoint(kind="tcp", pid=4242, tcp_port=12345)


@pytest.mark.parametrize(
    "content",
    [
        pytest.param("protocol=1\npid=4242\nunix=/tmp/a.sock\n", id="missing EOF"),
        pytest.param("protocol=1\npid=4242\nunix=/tmp/a.sockEOF\n", id="EOF not own line"),
        pytest.param(VALID_UNIX + "extra\n", id="content after EOF"),
        pytest.param("protocol=1\npid=4242\n", id="torn write"),
        pytest.param("EOF\n", id="empty body"),
        pytest.param("protocol=1\n\npid=4242\nunix=/tmp/a.sock\nEOF\n", id="empty line"),
        pytest.param("protocol=1\npid\nunix=/tmp/a.sock\nEOF\n", id="not key value"),
        pytest.param("protocol=1\npid=4242\nunix=/tmp/a.sock\nextra=1\nEOF\n", id="unknown key"),
        pytest.param("protocol=1\npid=4242\nmode=owner\nunix=/tmp/a.sock\nEOF\n", id="legacy mode key"),
        pytest.param("protocol=1\npid=4242\nunix=/tmp/a.sock\nauth_token=abc\nEOF\n", id="legacy auth key"),
        pytest.param("protocol=1\nprotocol=1\npid=4242\nunix=/tmp/a.sock\nEOF\n", id="duplicate key"),
        pytest.param("protocol=2\npid=4242\nunix=/tmp/a.sock\nEOF\n", id="unknown protocol"),
        pytest.param("protocol=abc\npid=4242\nunix=/tmp/a.sock\nEOF\n", id="non numeric protocol"),
        pytest.param("pid=4242\nunix=/tmp/a.sock\nEOF\n", id="missing protocol"),
        pytest.param("protocol=1\nunix=/tmp/a.sock\nEOF\n", id="missing pid"),
        pytest.param("protocol=1\npid=4242\nEOF\n", id="missing endpoint"),
        pytest.param("protocol=1\npid=4242\nunix=/tmp/a.sock\nsock=12345\nEOF\n", id="both endpoints"),
        pytest.param("protocol=1\npid=4242\nunix=tmp/a.sock\nEOF\n", id="relative unix path"),
        pytest.param(
            "protocol=1\npid=4242\nunix=/" + "a" * MAX_UNIX_PATH_BYTES + "\nEOF\n",
            id="oversize unix path",
        ),
        pytest.param("protocol=1\npid=4242\nsock=0\nEOF\n", id="port zero"),
        pytest.param("protocol=1\npid=4242\nsock=65536\nEOF\n", id="port range"),
        pytest.param("protocol=1\npid=4242\nsock=01234\nEOF\n", id="port leading zero"),
        pytest.param("protocol=1\npid=4242\nsock=12a45\nEOF\n", id="port not digits"),
        pytest.param("protocol=1\npid=0\nunix=/tmp/a.sock\nEOF\n", id="pid zero"),
        pytest.param("protocol=1\npid=-1\nunix=/tmp/a.sock\nEOF\n", id="pid negative"),
        pytest.param("protocol=1\npid=04242\nunix=/tmp/a.sock\nEOF\n", id="pid leading zero"),
        pytest.param("protocol=1\npid=42a42\nunix=/tmp/a.sock\nEOF\n", id="pid not digits"),
        pytest.param("protocol=1\npid=2147483648\nunix=/tmp/a.sock\nEOF\n", id="pid out of range"),
    ],
)
def test_parse_rejects_malformed(content):
    """TestCase: TestParseRejectsMalformed"""
    with pytest.raises(CoreGoError) as caught:
        parse_port_file(content.encode())

    assert caught.value.stage is Stage.PORT_FILE
    assert caught.value.reason is ReasonCode.PORT_FILE_INVALID


def test_parse_rejects_oversize_content():
    """TestCase: TestParseRejectsOversizeContent"""
    content = "protocol=1\npid=4242\nunix=/tmp/" + "a" * PORT_FILE_MAX_BYTES + "\nEOF\n"

    with pytest.raises(CoreGoError) as caught:
        parse_port_file(content.encode())

    assert caught.value.reason is ReasonCode.PORT_FILE_INVALID
    assert "exceeds" in caught.value.message


def test_parse_rejects_oversize_unix_path_exactly_at_limit():
    """unix 值恰好 103 字节可接受，104 字节拒绝。"""
    accepted = "/" + "a" * (MAX_UNIX_PATH_BYTES - 1)
    assert len(accepted.encode()) == MAX_UNIX_PATH_BYTES
    endpoint = parse_port_file(f"protocol=1\npid=1\nunix={accepted}\nEOF\n".encode())
    assert endpoint.unix_path == accepted

    rejected = "/" + "a" * MAX_UNIX_PATH_BYTES
    with pytest.raises(CoreGoError):
        parse_port_file(f"protocol=1\npid=1\nunix={rejected}\nEOF\n".encode())


def test_read_port_file_missing_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_port_file(tmp_path / "core.port")


def test_read_port_file_reads_valid_content(tmp_path):
    path = tmp_path / "core.port"
    path.write_bytes(VALID_UNIX.encode())

    assert read_port_file(path) == Endpoint(kind="unix", pid=4242, unix_path="/tmp/swanlab/core.sock")


def test_read_port_file_rejects_oversize_file(tmp_path):
    path = tmp_path / "core.port"
    path.write_bytes(b"a" * (PORT_FILE_MAX_BYTES + 1))

    with pytest.raises(CoreGoError) as caught:
        read_port_file(path)

    assert caught.value.reason is ReasonCode.PORT_FILE_INVALID


def test_endpoint_target():
    assert Endpoint(kind="unix", pid=1, unix_path="/tmp/core.sock").target == "unix:/tmp/core.sock"
    assert Endpoint(kind="tcp", pid=1, tcp_port=50051).target == "127.0.0.1:50051"
