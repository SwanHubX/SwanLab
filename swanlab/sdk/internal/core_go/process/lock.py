"""
@author: Nexisato
@file: lock.py
@time: 2026/11/16
@description: per-user runtime 目录的文件启动锁（core.lock）。

格式与 Go 侧 core/internal/pkg/portinfo/lock.go 保持一致：
  pid=<pid>
  created_at=<unix_ns>
  EOF

使用临时文件 + os.link 原子发布锁内容（不支持硬链接时回退 O_CREAT|O_EXCL），
支持死进程与超时陈旧锁回收；Release 仅删除归属本实例的锁文件。
"""

import ctypes
import os
import secrets
import stat
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple

from swanlab.sdk.internal.core_go.config import (
    LOCK_FILE_NAME,
    PORT_FILE_MAX_BYTES,
    START_LOCK_POLL_INTERVAL,
    START_LOCK_STALE_SECONDS,
    START_LOCK_TIMEOUT,
)
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.pkg import safe

_LOCK_FILE_MODE = 0o600
_NS_PER_SECOND = 1_000_000_000
_MAX_PID = (1 << 31) - 1
_EOF = "EOF"
_WINDOWS_QUERY_LIMITED = 0x1000
_WINDOWS_STILL_ACTIVE = 259
_WINDOWS_ERROR_ACCESS_DENIED = 5


def _lock_error(message: str, *, cause: Optional[BaseException] = None) -> CoreGoError:
    return CoreGoError(Stage.RUNTIME, ReasonCode.RUNTIME_UNAVAILABLE, message, cause=cause)


def is_process_alive(pid: int) -> bool:
    """跨平台判断指定 PID 的进程是否仍在运行。"""
    if pid <= 0:
        return False
    if os.name == "nt":
        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        handle = kernel32.OpenProcess(_WINDOWS_QUERY_LIMITED, False, int(pid))
        if not handle:
            return kernel32.GetLastError() == _WINDOWS_ERROR_ACCESS_DENIED
        try:
            exit_code = ctypes.c_ulong()
            if not kernel32.GetExitCodeProcess(handle, ctypes.byref(exit_code)):
                return False
            return exit_code.value == _WINDOWS_STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)

    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def _marshal_lock(pid: int, created_at: int) -> bytes:
    return f"pid={pid}\ncreated_at={created_at}\n{_EOF}\n".encode("ascii")


def _parse_lock(data: bytes) -> Tuple[int, int]:
    if len(data) > PORT_FILE_MAX_BYTES:
        raise ValueError(f"lock file exceeds {PORT_FILE_MAX_BYTES} bytes")
    text = data.decode("utf-8", errors="strict")
    content = text[:-1] if text.endswith("\n") else text
    suffix = f"\n{_EOF}"
    if not content.endswith(suffix):
        raise ValueError("lock file must end with an EOF line")
    body = content[: -len(suffix)]

    seen = {}
    for line in body.split("\n"):
        if not line:
            raise ValueError("lock file contains empty line")
        key, sep, value = line.partition("=")
        if not sep:
            raise ValueError("lock file line is not key=value")
        if key in seen:
            raise ValueError(f"lock file contains duplicate key {key!r}")
        if key not in ("pid", "created_at"):
            raise ValueError(f"lock file contains unknown key {key!r}")
        seen[key] = value

    if "pid" not in seen or "created_at" not in seen:
        raise ValueError("lock file missing required fields")

    pid_raw = seen["pid"]
    if not pid_raw.isdigit() or pid_raw.startswith("0"):
        raise ValueError(f"invalid pid {pid_raw!r}")
    pid = int(pid_raw)
    if pid < 1 or pid > _MAX_PID:
        raise ValueError(f"pid out of range: {pid}")

    ts_raw = seen["created_at"]
    if not ts_raw.isdigit():
        raise ValueError(f"invalid created_at {ts_raw!r}")
    created_at = int(ts_raw)
    if created_at <= 0:
        raise ValueError(f"created_at must be positive: {created_at}")

    return pid, created_at


def _read_lock_file(lock_path: Path) -> Tuple[int, int]:
    with lock_path.open("rb") as file:
        data = file.read(PORT_FILE_MAX_BYTES + 1)
    return _parse_lock(data)


def _create_lock_file(lock_path: Path, pid: int, created_at: int) -> bool:
    """原子创建锁文件；成功返回 True，已被占用返回 False。"""
    payload = _marshal_lock(pid, created_at)
    # Windows 下 os.open 默认文本模式会把 \n 转成 \r\n，必须显式附加 O_BINARY。
    open_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    tmp_path = lock_path.parent / f".{LOCK_FILE_NAME}-{pid}-{secrets.token_hex(6)}"
    try:
        fd = os.open(str(tmp_path), open_flags, _LOCK_FILE_MODE)
        try:
            view = memoryview(payload)
            written = 0
            while written < len(view):
                written += os.write(fd, view[written:])
        finally:
            os.close(fd)
        if os.name != "nt":
            os.chmod(tmp_path, _LOCK_FILE_MODE)
        try:
            os.link(str(tmp_path), str(lock_path))
            return True
        except FileExistsError:
            return False
        except OSError:
            # 少数不支持硬链接的文件系统回退到 O_CREAT|O_EXCL。
            pass
    finally:
        with safe.block(OSError, message=None):
            os.unlink(str(tmp_path))

    try:
        fd = os.open(str(lock_path), open_flags, _LOCK_FILE_MODE)
    except FileExistsError:
        return False

    try:
        view = memoryview(payload)
        written = 0
        while written < len(view):
            written += os.write(fd, view[written:])
    except BaseException:
        os.close(fd)
        with safe.block(OSError, message=None):
            os.unlink(str(lock_path))
        raise
    os.close(fd)
    return True


def _reclaim_if_stale(
    lock_path: Path,
    now_ns: int,
    stale_ns: int,
    alive_fn: Callable[[int], bool],
) -> bool:
    try:
        info = lock_path.lstat()
    except FileNotFoundError:
        return True
    except OSError as exc:
        raise _lock_error(f"cannot stat lock file {lock_path}: {exc}", cause=exc) from exc

    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise _lock_error(f"lock path {lock_path} is not a regular file")

    try:
        holder_pid, created_at = _read_lock_file(lock_path)
    except FileNotFoundError:
        return True
    except (OSError, ValueError):
        # 尚未写完的锁文件仅在超过陈旧阈值后才按崩溃残留清理。
        if (now_ns - int(info.st_mtime * _NS_PER_SECOND)) < stale_ns:
            return False
    else:
        expired = (now_ns - created_at) >= stale_ns
        if not expired and alive_fn(holder_pid):
            return False

    try:
        latest = lock_path.lstat()
    except FileNotFoundError:
        return True
    except OSError as exc:
        raise _lock_error(f"cannot re-stat lock file {lock_path}: {exc}", cause=exc) from exc

    if (latest.st_dev, latest.st_ino) != (info.st_dev, info.st_ino):
        return False

    try:
        os.unlink(str(lock_path))
    except FileNotFoundError:
        return True
    except OSError as exc:
        raise _lock_error(f"cannot remove stale lock {lock_path}: {exc}", cause=exc) from exc
    return True


@dataclass
class StartupLock:
    """已持有的启动锁句柄，支持幂等 release 与 with 上下文。"""

    path: Path
    pid: int
    created_at: int
    _released: bool = False

    def release(self) -> None:
        if self._released:
            return
        self._released = True

        try:
            info = self.path.lstat()
        except FileNotFoundError:
            return
        except OSError:
            return
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
            return

        try:
            current_pid, current_created = _read_lock_file(self.path)
        except (OSError, ValueError):
            return
        if current_pid != self.pid or current_created != self.created_at:
            return

        with safe.block(OSError, message=None):
            os.unlink(str(self.path))

    def __enter__(self) -> "StartupLock":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.release()


def try_acquire_start_lock(
    runtime_dir: Path,
    *,
    validate_dir: Optional[Callable[[Path], None]] = None,
    stale_seconds: float = START_LOCK_STALE_SECONDS,
    now_ns_fn: Optional[Callable[[], int]] = None,
    is_alive_fn: Optional[Callable[[int], bool]] = None,
) -> Optional[StartupLock]:
    """尝试非阻塞获取启动锁；若被存活持有者占用返回 None。"""
    if validate_dir is not None:
        validate_dir(runtime_dir)

    lock_path = Path(runtime_dir) / LOCK_FILE_NAME
    clock = now_ns_fn or time.time_ns
    alive_fn = is_alive_fn or is_process_alive
    stale_ns = max(0, int(stale_seconds * _NS_PER_SECOND))
    pid = os.getpid()
    created_at = clock()

    try:
        if _create_lock_file(lock_path, pid, created_at):
            return StartupLock(path=lock_path, pid=pid, created_at=created_at)
    except OSError as exc:
        raise _lock_error(f"failed to create startup lock {lock_path}: {exc}", cause=exc) from exc

    if not _reclaim_if_stale(lock_path, created_at, stale_ns, alive_fn):
        return None

    created_at = clock()
    try:
        if _create_lock_file(lock_path, pid, created_at):
            return StartupLock(path=lock_path, pid=pid, created_at=created_at)
    except OSError as exc:
        raise _lock_error(f"failed to create startup lock {lock_path}: {exc}", cause=exc) from exc
    return None


def acquire_start_lock(
    runtime_dir: Path,
    *,
    validate_dir: Optional[Callable[[Path], None]] = None,
    timeout: float = START_LOCK_TIMEOUT,
    stale_seconds: float = START_LOCK_STALE_SECONDS,
    poll_interval: float = START_LOCK_POLL_INTERVAL,
    now_ns_fn: Optional[Callable[[], int]] = None,
    is_alive_fn: Optional[Callable[[int], bool]] = None,
    sleep_fn: Optional[Callable[[float], None]] = None,
) -> StartupLock:
    """在 runtime 目录竞争启动锁，直到获取成功或超时抛出 CoreGoError。"""
    clock = now_ns_fn or time.time_ns
    sleeper = sleep_fn or time.sleep
    deadline_ns = clock() + max(0, int(timeout * _NS_PER_SECOND))

    while True:
        lock = try_acquire_start_lock(
            runtime_dir,
            validate_dir=validate_dir,
            stale_seconds=stale_seconds,
            now_ns_fn=clock,
            is_alive_fn=is_alive_fn,
        )
        if lock is not None:
            return lock

        remaining_ns = deadline_ns - clock()
        if remaining_ns <= 0:
            raise _lock_error(f"timed out waiting for core startup lock: {Path(runtime_dir) / LOCK_FILE_NAME}")
        sleeper(min(poll_interval, remaining_ns / _NS_PER_SECOND))
