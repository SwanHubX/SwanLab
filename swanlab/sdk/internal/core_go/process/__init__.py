"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: 子进程与 runtime 租约。

负责私有 runtime 目录、Popen 启动、等待端点、进程早死检测与
wait → terminate → kill 的分级回收。不构造 protobuf，不解释 session_id。

约束：
- runtime 固定放在系统临时根下，不使用 cwd / TMPDIR / log_dir；
- spawn 只负责拉起进程，端点归属由 discovery 校验；
- 回收与清理永不抛出，避免覆盖调用方的原始异常。
"""

import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.core_go.config import (
    ENDPOINT_POLL_INTERVAL,
    KILL_TIMEOUT,
    MAX_UNIX_PATH_BYTES,
    PORT_FILE_NAME,
    PORT_FILE_TIMEOUT,
    RUNTIME_DIR_MODE,
    RUNTIME_DIR_PREFIX,
    STDERR_FILE_NAME,
    STDERR_TAIL_BYTES,
    TERMINATE_TIMEOUT,
    WAIT_TIMEOUT,
)
from swanlab.sdk.internal.core_go.discovery import Endpoint, read_port_file, validate_ownership
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.pkg import console, fork

__all__ = [
    "RuntimeLease",
    "acquire_runtime",
    "cleanup_runtime",
    "validate_core_binary",
    "spawn_core",
    "wait_endpoint",
    "reap_process",
    "stderr_tail",
]

_SOCKET_NAME = "core.sock"
_EXIT_ENDPOINT_BUSY = 3
_EXIT_USAGE_ERROR = 2


@dataclass
class RuntimeLease:
    """一次 core 进程独占的 runtime 目录。

    owner_pid 是创建租约的 Python 进程；fork child 据此拒绝回收父进程的资源。
    """

    dir: Path
    port_file: Path
    stderr_file: Path
    owner_pid: int


def _runtime_root() -> Path:
    """返回 runtime 根目录：POSIX 固定 /tmp，Windows 使用 %TEMP%。"""
    if os.name == "nt":
        return Path(os.environ.get("TEMP") or os.environ.get("TMP") or tempfile.gettempdir())
    return Path("/tmp")


def _runtime_error(message: str, *, detail: Optional[str] = None, cause: Optional[BaseException] = None):
    return CoreGoError(Stage.RUNTIME, ReasonCode.RUNTIME_UNAVAILABLE, message, detail=detail, cause=cause)


def _validate_runtime_dir(path: Path) -> None:
    """校验 runtime 是当前用户拥有的 0700 真实目录，拒绝符号链接。"""
    try:
        info = path.lstat()
    except OSError as exc:
        raise _runtime_error(f"cannot stat runtime dir {path}: {exc}", cause=exc) from exc

    if stat.S_ISLNK(info.st_mode):
        raise _runtime_error(f"runtime dir {path} is a symlink")
    if not stat.S_ISDIR(info.st_mode):
        raise _runtime_error(f"runtime path {path} is not a directory")
    # Windows 的 chmod 只切换只读位，没有 POSIX 0700/uid 语义。
    if os.name != "nt":
        if stat.S_IMODE(info.st_mode) != RUNTIME_DIR_MODE:
            raise _runtime_error(f"runtime dir {path} mode is {oct(stat.S_IMODE(info.st_mode))}, want 0o700")
        if hasattr(os, "getuid") and info.st_uid != os.getuid():
            raise _runtime_error(f"runtime dir {path} is not owned by current user")


def acquire_runtime() -> RuntimeLease:
    """创建一次性的私有 runtime 目录。

    目录名带固定前缀，权限 0700；UDS 路径超长时直接失败，不回退 cwd。
    """
    root = _runtime_root()
    try:
        runtime_dir = Path(tempfile.mkdtemp(prefix=RUNTIME_DIR_PREFIX, dir=str(root)))
    except OSError as exc:
        raise _runtime_error(f"cannot create runtime dir under {root}: {exc}", cause=exc) from exc

    try:
        os.chmod(runtime_dir, RUNTIME_DIR_MODE)
        _validate_runtime_dir(runtime_dir)
        if os.name != "nt":
            socket_path = runtime_dir / _SOCKET_NAME
            if len(os.fsencode(str(socket_path))) > MAX_UNIX_PATH_BYTES:
                raise _runtime_error(f"unix socket path exceeds {MAX_UNIX_PATH_BYTES} bytes: {socket_path}")
    except BaseException:
        shutil.rmtree(runtime_dir, ignore_errors=True)
        raise

    return RuntimeLease(
        dir=runtime_dir,
        port_file=runtime_dir / PORT_FILE_NAME,
        stderr_file=runtime_dir / STDERR_FILE_NAME,
        owner_pid=fork.current_pid(),
    )


def cleanup_runtime(lease: RuntimeLease) -> None:
    """删除本进程创建的 runtime 目录；任何失败只记录，不抛出。"""
    try:
        if fork.is_forked(lease.owner_pid):
            return
        runtime_dir = Path(lease.dir)
        if not runtime_dir.name.startswith(RUNTIME_DIR_PREFIX):
            return
        if runtime_dir.is_symlink() or not runtime_dir.is_dir():
            return
        if hasattr(os, "getuid") and runtime_dir.lstat().st_uid != os.getuid():
            return
        shutil.rmtree(runtime_dir)
    except Exception as exc:  # noqa: BLE001 - 清理失败不能覆盖原始异常
        console.warning(f"failed to cleanup core runtime {lease.dir}: {exc}")


def stderr_tail(lease: RuntimeLease) -> str:
    """读取有界的 stderr 尾部，文件不存在时返回空串。"""
    try:
        with Path(lease.stderr_file).open("rb") as file:
            file.seek(0, os.SEEK_END)
            size = file.tell()
            file.seek(max(0, size - STDERR_TAIL_BYTES))
            data = file.read(STDERR_TAIL_BYTES)
    except OSError:
        return ""
    return data.decode("utf-8", "replace")


def validate_core_binary(binary: Path) -> None:
    """校验 core 二进制存在且可执行；launcher 在创建 runtime 前复用它。"""
    binary_path = Path(binary)
    if not binary_path.is_file():
        raise CoreGoError(
            Stage.SPAWN,
            ReasonCode.BINARY_NOT_FOUND,
            f"core binary not found: {binary_path}",
        )
    if not os.access(binary_path, os.X_OK):
        raise CoreGoError(
            Stage.SPAWN,
            ReasonCode.SPAWN_FAILED,
            f"core binary is not executable: {binary_path}",
        )


def spawn_core(binary: Path, lease: RuntimeLease, parent_pid: Optional[int] = None) -> subprocess.Popen:
    """拉起 swanlab-core，返回 Popen。

    只传 --port-filename / --parent-pid；不传 --listen / --detach，
    避免外部环境变量改变本次启动行为。子进程 cwd 是 runtime 目录，
    因此相对路径必须先按调用方 cwd 解析为绝对路径。
    """
    try:
        binary_path = Path(binary).resolve()
    except OSError as exc:
        raise CoreGoError(
            Stage.SPAWN,
            ReasonCode.SPAWN_FAILED,
            f"failed to resolve core binary {binary}: {exc}",
            cause=exc,
        ) from exc
    validate_core_binary(binary_path)

    env = dict(os.environ)
    env.pop("SWANLAB_CORE_LISTEN", None)
    env.pop("SWANLAB_CORE_PARENT_PID", None)

    argv = [
        str(binary_path),
        "--port-filename",
        str(lease.port_file),
        "--parent-pid",
        str(parent_pid if parent_pid is not None else fork.current_pid()),
    ]

    kwargs = {
        "args": argv,
        "cwd": str(lease.dir),
        "env": env,
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
    }
    if sys.platform == "win32":
        kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW
    else:
        kwargs["start_new_session"] = True

    try:
        with Path(lease.stderr_file).open("ab") as stderr:
            return subprocess.Popen(stderr=stderr, **kwargs)
    except OSError as exc:
        raise CoreGoError(
            Stage.SPAWN,
            ReasonCode.SPAWN_FAILED,
            f"failed to spawn core: {exc}",
            detail=stderr_tail(lease),
            cause=exc,
        ) from exc


def _early_exit_error(proc: subprocess.Popen, lease: RuntimeLease, exit_code: Optional[int]) -> CoreGoError:
    if exit_code == _EXIT_ENDPOINT_BUSY:
        reason = ReasonCode.ENDPOINT_BUSY
        message = "core exited because the endpoint is already in use"
    elif exit_code == _EXIT_USAGE_ERROR:
        reason = ReasonCode.USAGE_ERROR
        message = "core exited with usage error"
    else:
        reason = ReasonCode.EARLY_EXIT
        message = f"core exited before publishing endpoint (exit code {exit_code})"
    return CoreGoError(Stage.PORT_FILE, reason, message, detail=stderr_tail(lease))


def wait_endpoint(proc: subprocess.Popen, lease: RuntimeLease, timeout: float = PORT_FILE_TIMEOUT) -> Endpoint:
    """轮询 port-file 直到解析成功，并校验端点归属。

    子进程早死时按退出码分类；解析成功但 pid 不匹配立即失败，
    因为 Go 不会重写该文件。
    """
    deadline = time.monotonic() + timeout
    while True:
        exit_code = proc.poll()
        if exit_code is not None:
            raise _early_exit_error(proc, lease, exit_code)

        try:
            endpoint = read_port_file(lease.port_file)
        except FileNotFoundError:
            endpoint = None
        except CoreGoError as exc:
            if exc.reason is not ReasonCode.PORT_FILE_INVALID:
                raise
            endpoint = None

        if endpoint is not None:
            validate_ownership(endpoint, proc.pid, lease.dir)
            return endpoint

        if time.monotonic() >= deadline:
            exit_code = proc.poll()
            if exit_code is not None:
                raise _early_exit_error(proc, lease, exit_code)
            raise CoreGoError(
                Stage.PORT_FILE,
                ReasonCode.PORT_FILE_TIMEOUT,
                f"port-file not ready within {timeout:.1f}s",
                detail=stderr_tail(lease),
            )
        time.sleep(ENDPOINT_POLL_INTERVAL)


def reap_process(proc: Optional[subprocess.Popen]) -> Optional[int]:
    """分级回收子进程：wait(10) → terminate + wait(3) → kill + wait(2)。

    返回最后一次 poll() 的退出码；所有 ProcessLookupError / TimeoutExpired
    都被吞掉，函数本身永不抛出。
    """
    if proc is None:
        return None
    try:
        return proc.wait(timeout=WAIT_TIMEOUT)
    except subprocess.TimeoutExpired:
        pass
    except OSError:
        return proc.poll()

    try:
        proc.terminate()
        return proc.wait(timeout=TERMINATE_TIMEOUT)
    except subprocess.TimeoutExpired:
        pass
    except OSError:
        pass

    try:
        proc.kill()
        return proc.wait(timeout=KILL_TIMEOUT)
    except (subprocess.TimeoutExpired, OSError):
        return proc.poll()
