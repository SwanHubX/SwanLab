"""
@author: Nexisato
@file: portfile.py
@time: 2026/10/27 15:40
@description: port-file 解析、端点发现与归属校验。

解析 port-file v1，产出 Endpoint（unix / tcp 二选一）并校验 pid 与目录归属。
纯函数，不建立 channel，也不启动进程。

本地协议的唯一事实源是 core/internal/pkg/portinfo/portinfo.go：
本模块逐条复刻 Parse 的接受与拒绝规则，不额外放宽。

约束：
- 只做纯解析与校验，不 spawn、不写文件、不 import grpc；
- 文件不存在时抛 FileNotFoundError，交给轮询方决定是否继续等待；
- 解析失败统一抛 CoreGoError(stage=PORT_FILE, reason=PORT_FILE_INVALID)。
"""

import socket
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from swanlab.sdk.internal.core_go.config import (
    ENDPOINT_PROBE_TIMEOUT,
    MAX_UNIX_PATH_BYTES,
    PORT_FILE_EOF,
    PORT_FILE_MAX_BYTES,
    CoreGoError,
    ReasonCode,
    Stage,
)

__all__ = [
    "Endpoint",
    "discover_endpoint",
    "parse_port_file",
    "probe_endpoint",
    "read_port_file",
    "validate_discovered_endpoint",
    "validate_ownership",
]

_MAX_PID = 1 << 31 - 1
_MAX_PORT = 65535
_LINE_ECHO_LIMIT = 120


@dataclass(frozen=True)
class Endpoint:
    """port-file 描述的服务端点。

    kind 为 ``"unix"`` 时使用 ``unix_path``，为 ``"tcp"`` 时使用 ``tcp_port``。
    不保存 port-file 原始文本，避免被误当成长期注册表。
    """

    kind: str
    pid: int
    unix_path: Optional[str] = None
    tcp_port: Optional[int] = None

    @property
    def target(self) -> str:
        """gRPC target 字符串：``unix:/abs/path`` 或 ``127.0.0.1:port``。"""
        if self.kind == "unix":
            return f"unix:{self.unix_path}"
        return f"127.0.0.1:{self.tcp_port}"


def _invalid(message: str) -> CoreGoError:
    return CoreGoError(Stage.PORT_FILE, ReasonCode.PORT_FILE_INVALID, message)


def _echo(raw: bytes) -> str:
    """只回显被拒绝的内容，并截断，避免大文件进入日志。"""
    return repr(raw.decode("utf-8", "replace")[:_LINE_ECHO_LIMIT])


def _parse_decimal(value: bytes, field: str, upper: int) -> int:
    """解析十进制字段：非空、全数字、无前导零、范围 1..upper。"""
    if not value or value.startswith(b"0") or not value.isdigit():
        raise _invalid(f"invalid {field} {_echo(value)}")
    number = int(value)
    if number < 1 or number > upper:
        raise _invalid(f"invalid {field} {_echo(value)}")
    return number


def parse_port_file(data: bytes) -> Endpoint:
    """解析 port-file v1 内容。

    与 Go portinfo.Parse 保持一致：键值行顺序不限，独立 EOF 行结尾，
    末尾最多一个换行；重复键、未知键、空行、缺失字段都拒绝。
    """
    if len(data) > PORT_FILE_MAX_BYTES:
        raise _invalid(f"port-file exceeds {PORT_FILE_MAX_BYTES} bytes")

    content = data[:-1] if data.endswith(b"\n") else data
    if content == PORT_FILE_EOF:
        body = b""
    elif content.endswith(b"\n" + PORT_FILE_EOF):
        body = content[: -(len(PORT_FILE_EOF) + 1)]
    else:
        raise _invalid("port-file must end with an EOF line")

    seen: Dict[bytes, bool] = {}
    values: Dict[bytes, bytes] = {}
    for line in body.split(b"\n"):
        if not line:
            raise _invalid("port-file contains empty line")
        key, sep, value = line.partition(b"=")
        if not sep:
            raise _invalid(f"port-file line is not key=value: {_echo(line)}")
        if key in seen:
            raise _invalid(f"port-file contains duplicate key {_echo(key)}")
        if key not in (b"protocol", b"pid", b"unix", b"sock"):
            raise _invalid(f"port-file contains unknown key {_echo(key)}")
        seen[key] = True
        values[key] = value

    if b"protocol" not in seen:
        raise _invalid('port-file missing "protocol" field')
    if values[b"protocol"] != b"1":
        raise _invalid(f"unsupported port-file protocol {_echo(values[b'protocol'])}")
    if b"pid" not in seen:
        raise _invalid('port-file missing "pid" field')

    pid = _parse_decimal(values[b"pid"], "pid", _MAX_PID)

    has_unix = b"unix" in seen
    has_sock = b"sock" in seen
    if has_unix == has_sock:
        raise _invalid('port-file must contain exactly one of "unix" or "sock"')

    if has_unix:
        unix_bytes = values[b"unix"]
        if not unix_bytes.startswith(b"/") or len(unix_bytes) > MAX_UNIX_PATH_BYTES:
            raise _invalid("unix endpoint must be an absolute path within length limit")
        return Endpoint(kind="unix", pid=pid, unix_path=unix_bytes.decode("utf-8", "surrogateescape"))

    port = _parse_decimal(values[b"sock"], "sock port", _MAX_PORT)
    return Endpoint(kind="tcp", pid=pid, tcp_port=port)


def read_port_file(path: Path) -> Endpoint:
    """读取并解析 port-file。

    最多读 ``PORT_FILE_MAX_BYTES + 1`` 字节，用多出的 1 字节判超长。
    文件不存在时保留 FileNotFoundError，表示「还没写好」。
    """
    with Path(path).open("rb") as file:
        data = file.read(PORT_FILE_MAX_BYTES + 1)
    return parse_port_file(data)


def validate_discovered_endpoint(endpoint: Endpoint, runtime_dir: Path) -> None:
    """校验待复用端点落在指定的 runtime 目录边界内（不绑定创建者 PID）。"""
    if endpoint.kind == "unix":
        socket_path = Path(endpoint.unix_path or "")
        if socket_path.parent != Path(runtime_dir):
            raise CoreGoError(
                Stage.PORT_FILE,
                ReasonCode.OWNERSHIP_MISMATCH,
                f"unix socket {socket_path} is outside runtime dir {runtime_dir}",
            )
        if socket_path.is_symlink():
            raise CoreGoError(
                Stage.PORT_FILE,
                ReasonCode.OWNERSHIP_MISMATCH,
                f"unix socket {socket_path} is a symlink",
            )
    elif endpoint.kind == "tcp":
        if endpoint.tcp_port is None or not 1 <= endpoint.tcp_port <= _MAX_PORT:
            raise CoreGoError(
                Stage.PORT_FILE,
                ReasonCode.OWNERSHIP_MISMATCH,
                f"invalid tcp port {endpoint.tcp_port}",
            )
    else:
        raise CoreGoError(
            Stage.PORT_FILE,
            ReasonCode.OWNERSHIP_MISMATCH,
            f"unsupported endpoint kind {endpoint.kind!r}",
        )


def validate_ownership(endpoint: Endpoint, process_pid: int, runtime_dir: Path) -> None:
    """校验端点确实由本次 spawn 的进程创建。

    - pid 必须等于 Popen 的 pid；
    - UDS 必须落在本次 runtime 目录内，且 socket 本身不是符号链接；
    - TCP 只允许 loopback，port-file 没有 host 字段，无需额外解析。
    """
    if endpoint.pid != process_pid:
        raise CoreGoError(
            Stage.PORT_FILE,
            ReasonCode.OWNERSHIP_MISMATCH,
            f"port-file pid {endpoint.pid} does not match spawned pid {process_pid}",
        )
    validate_discovered_endpoint(endpoint, runtime_dir)


def probe_endpoint(endpoint: Endpoint, timeout: float = ENDPOINT_PROBE_TIMEOUT) -> bool:
    """用短超时探测端点是否仍有 core 进程在监听。"""
    if endpoint.kind == "unix":
        if not endpoint.unix_path or not hasattr(socket, "AF_UNIX"):
            return False
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            sock.settimeout(timeout)
            sock.connect(endpoint.unix_path)
            return True
        except OSError:
            return False
        finally:
            sock.close()

    if endpoint.kind == "tcp":
        if endpoint.tcp_port is None or not 1 <= endpoint.tcp_port <= _MAX_PORT:
            return False
        try:
            with socket.create_connection(("127.0.0.1", endpoint.tcp_port), timeout=timeout):
                return True
        except OSError:
            return False

    return False


def discover_endpoint(
    port_file: Path,
    runtime_dir: Path,
    *,
    probe_timeout: float = ENDPOINT_PROBE_TIMEOUT,
) -> Optional[Endpoint]:
    """发现并探测已有可用端点；缺失、损坏、越界或不可达时返回 None。"""
    try:
        endpoint = read_port_file(port_file)
    except OSError:
        return None
    except CoreGoError as exc:
        if exc.reason is ReasonCode.PORT_FILE_INVALID:
            return None
        raise

    try:
        validate_discovered_endpoint(endpoint, runtime_dir)
    except CoreGoError:
        return None

    if not probe_endpoint(endpoint, timeout=probe_timeout):
        return None
    return endpoint
