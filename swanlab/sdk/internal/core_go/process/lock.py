"""
@author: Nexisato
@file: lock.py
@time: 2026/11/16
@description: per-user runtime 目录的内核独占启动锁（core.lock）。

POSIX 使用 fcntl.flock(LOCK_EX | LOCK_NB)，Windows 使用 LockFileEx，
与 Go 侧 core/internal/pkg/portinfo/lock.go 共享同一把内核锁：
- 不删除 core.lock 文件本身，从根本上消除检查→删除→重建的 ABA 竞态；
- 持有者进程崩溃或被 kill -9 后，内核自动释放锁；
- 持锁期间写入 pid / created_at 元数据供诊断查看。
"""

import ctypes
import errno
import os
import stat
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

from swanlab.sdk.internal.core_go.config import (
    LOCK_FILE_NAME,
    START_LOCK_POLL_INTERVAL,
    START_LOCK_TIMEOUT,
)
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage
from swanlab.sdk.internal.pkg import safe

if os.name != "nt":
    import fcntl as _fcntl
else:  # pragma: no cover - Windows 分支走 LockFileEx
    _fcntl: Any = None

_LOCK_FILE_MODE = 0o600
_NS_PER_SECOND = 1_000_000_000
_EOF = "EOF"
_LOCKFILE_FAIL_IMMEDIATELY = 0x00000001
_LOCKFILE_EXCLUSIVE_LOCK = 0x00000002
_ERROR_LOCK_VIOLATION = 33
_LOCK_RANGE_BYTES = 1


def _lock_error(message: str, *, cause: Optional[BaseException] = None) -> CoreGoError:
    return CoreGoError(Stage.RUNTIME, ReasonCode.RUNTIME_UNAVAILABLE, message, cause=cause)


def _marshal_lock(pid: int, created_at: int) -> bytes:
    return f"pid={pid}\ncreated_at={created_at}\n{_EOF}\n".encode("ascii")


def _write_all(fd: int, data: bytes) -> None:
    view = memoryview(data)
    written = 0
    while written < len(view):
        written += os.write(fd, view[written:])


def _write_lock_metadata(fd: int, pid: int, created_at: int) -> None:
    os.lseek(fd, 0, os.SEEK_SET)
    os.ftruncate(fd, 0)
    _write_all(fd, _marshal_lock(pid, created_at))


if os.name == "nt":  # pragma: no cover - Windows 内核文件锁
    import msvcrt

    class _OVERLAPPED(ctypes.Structure):
        _fields_ = [
            ("Internal", ctypes.c_void_p),
            ("InternalHigh", ctypes.c_void_p),
            ("Offset", ctypes.c_uint32),
            ("OffsetHigh", ctypes.c_uint32),
            ("hEvent", ctypes.c_void_p),
        ]

    _kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)  # type: ignore[attr-defined]
    _kernel32.LockFileEx.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.POINTER(_OVERLAPPED),
    ]
    _kernel32.LockFileEx.restype = ctypes.c_int
    _kernel32.UnlockFileEx.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.POINTER(_OVERLAPPED),
    ]
    _kernel32.UnlockFileEx.restype = ctypes.c_int

    def _try_flock_exclusive(fd: int) -> bool:
        handle = msvcrt.get_osfhandle(fd)
        ov = _OVERLAPPED()
        ok = _kernel32.LockFileEx(
            handle,
            _LOCKFILE_EXCLUSIVE_LOCK | _LOCKFILE_FAIL_IMMEDIATELY,
            0,
            _LOCK_RANGE_BYTES,
            0,
            ctypes.byref(ov),
        )
        if ok:
            return True
        err = ctypes.get_last_error()
        if err == _ERROR_LOCK_VIOLATION:
            return False
        raise OSError(err, ctypes.FormatError(err))

    def _unlock_fd(fd: int) -> None:
        try:
            handle = msvcrt.get_osfhandle(fd)
            ov = _OVERLAPPED()
            _kernel32.UnlockFileEx(handle, 0, _LOCK_RANGE_BYTES, 0, ctypes.byref(ov))
        except OSError:
            pass
        finally:
            with safe.block(OSError, message=None):
                os.close(fd)

else:

    def _try_flock_exclusive(fd: int) -> bool:
        try:
            _fcntl.flock(fd, _fcntl.LOCK_EX | _fcntl.LOCK_NB)
            return True
        except BlockingIOError:
            return False
        except OSError as exc:
            if exc.errno in (errno.EAGAIN, errno.EWOULDBLOCK, errno.EACCES):
                return False
            raise

    def _unlock_fd(fd: int) -> None:
        with safe.block(OSError, message=None):
            _fcntl.flock(fd, _fcntl.LOCK_UN)
        with safe.block(OSError, message=None):
            os.close(fd)


def _open_lock_fd(lock_path: Path) -> int:
    if os.name == "nt":
        try:
            info = lock_path.lstat()
        except FileNotFoundError:
            pass
        except OSError as exc:
            raise _lock_error(f"cannot stat lock file {lock_path}: {exc}", cause=exc) from exc
        else:
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
                raise _lock_error(f"lock path {lock_path} is not a regular file")

    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(str(lock_path), flags, _LOCK_FILE_MODE)
    except OSError as exc:
        raise _lock_error(f"failed to open startup lock {lock_path}: {exc}", cause=exc) from exc

    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            raise _lock_error(f"lock path {lock_path} is not a regular file")
    except BaseException:
        with safe.block(OSError, message=None):
            os.close(fd)
        raise
    return fd


@dataclass
class StartupLock:
    """已持有的内核启动锁句柄，支持幂等 release 与 with 上下文。"""

    path: Path
    pid: int
    created_at: int
    _fd: Optional[int] = None
    _mu: threading.Lock = field(default_factory=threading.Lock, repr=False, compare=False)

    def release(self) -> None:
        with self._mu:
            fd = self._fd
            if fd is None:
                return
            self._fd = None
        _unlock_fd(fd)

    def __enter__(self) -> "StartupLock":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.release()


def try_acquire_start_lock(
    runtime_dir: Path,
    *,
    validate_dir: Optional[Callable[[Path], None]] = None,
    now_ns_fn: Optional[Callable[[], int]] = None,
) -> Optional[StartupLock]:
    """尝试非阻塞获取内核独占启动锁；若被占用返回 None。"""
    if validate_dir is not None:
        validate_dir(runtime_dir)

    lock_path = Path(runtime_dir) / LOCK_FILE_NAME
    fd = _open_lock_fd(lock_path)

    try:
        if not _try_flock_exclusive(fd):
            with safe.block(OSError, message=None):
                os.close(fd)
            return None

        clock = now_ns_fn or time.time_ns
        pid = os.getpid()
        created_at = clock()
        _write_lock_metadata(fd, pid, created_at)
        return StartupLock(path=lock_path, pid=pid, created_at=created_at, _fd=fd)
    except CoreGoError:
        _unlock_fd(fd)
        raise
    except OSError as exc:
        _unlock_fd(fd)
        raise _lock_error(f"failed to lock {lock_path}: {exc}", cause=exc) from exc
    except BaseException:
        _unlock_fd(fd)
        raise


def acquire_start_lock(
    runtime_dir: Path,
    *,
    validate_dir: Optional[Callable[[Path], None]] = None,
    timeout: float = START_LOCK_TIMEOUT,
    poll_interval: float = START_LOCK_POLL_INTERVAL,
    now_ns_fn: Optional[Callable[[], int]] = None,
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
            now_ns_fn=clock,
        )
        if lock is not None:
            return lock

        remaining_ns = deadline_ns - clock()
        if remaining_ns <= 0:
            raise _lock_error(f"timed out waiting for core startup lock: {Path(runtime_dir) / LOCK_FILE_NAME}")
        sleeper(min(poll_interval, remaining_ns / _NS_PER_SECOND))
