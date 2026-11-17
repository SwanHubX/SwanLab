"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: 本地运行时环境：per-user 目录租约、二进制定位、内核启动锁与端点发现。

不负责启动子进程（Popen）也不导入 gRPC；只处理本地文件系统契约与端点探测。
"""

import os
import shutil
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Optional

from swanlab.sdk.internal.core_go.config import (
    BINARY_NAME,
    CORE_BINARY_ENV,
    MAX_UNIX_PATH_BYTES,
    PORT_FILE_NAME,
    RUNTIME_DIR_MODE,
    RUNTIME_DIR_PREFIX,
    SOCKET_NAME,
    START_LOCK_POLL_INTERVAL,
    START_LOCK_TIMEOUT,
    STDERR_FILE_NAME,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.runtime.lock import StartupLock
from swanlab.sdk.internal.core_go.runtime.lock import (
    acquire_start_lock as _acquire_start_lock,
)
from swanlab.sdk.internal.core_go.runtime.lock import (
    try_acquire_start_lock as _try_acquire_start_lock,
)
from swanlab.sdk.internal.core_go.runtime.portfile import (
    Endpoint,
    discover_endpoint,
    parse_port_file,
    probe_endpoint,
    read_port_file,
    validate_discovered_endpoint,
    validate_ownership,
)
from swanlab.sdk.internal.pkg import console, fork

__all__ = [
    "Endpoint",
    "RuntimeLease",
    "StartupLock",
    "acquire_runtime",
    "acquire_start_lock",
    "cleanup_runtime",
    "default_runtime_dir",
    "discover_endpoint",
    "locate_core_binary",
    "parse_port_file",
    "probe_endpoint",
    "read_port_file",
    "try_acquire_start_lock",
    "validate_core_binary",
    "validate_discovered_endpoint",
    "validate_ownership",
]

_WINDOWS_USER_ENV_KEYS = ("USERNAME", "USER", "LOGNAME", "LNAME")
# swanlab/sdk/internal/core_go/runtime/__init__.py -> swanlab/bin/
_PACKAGED_BINARY = (
    Path(__file__).resolve().parents[4] / "bin" / (f"{BINARY_NAME}.exe" if os.name == "nt" else BINARY_NAME)
)


@dataclass
class RuntimeLease:
    """core 进程的 runtime 目录租约。

    owner_pid 是创建租约的 Python 进程；fork child 据此拒绝回收父进程的资源。
    """

    dir: Path
    port_file: Path
    stderr_file: Path
    owner_pid: int


def locate_core_binary() -> Path:
    """定位 swanlab-core：SWANLAB_CORE_BINARY 优先，其次包内 swanlab/bin。"""
    override = os.environ.get(CORE_BINARY_ENV, "").strip()
    binary = Path(override).expanduser() if override else _PACKAGED_BINARY
    if not binary.is_file():
        raise CoreGoError(
            Stage.SPAWN,
            ReasonCode.BINARY_NOT_FOUND,
            f"swanlab-core binary not found: {binary}; "
            f"set {CORE_BINARY_ENV} or install a platform wheel that bundles swanlab-core",
        )
    return binary.resolve()


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


def _runtime_root() -> Path:
    """返回 runtime 根目录：POSIX 固定 /tmp，Windows 使用 %TEMP%。"""
    if os.name == "nt":
        return Path(os.environ.get("TEMP", "").strip() or os.environ.get("TMP", "").strip() or tempfile.gettempdir())
    return Path("/tmp")


def _sanitize_user_segment(raw: str) -> str:
    cleaned = "".join(
        ch if (("a" <= ch <= "z") or ("A" <= ch <= "Z") or ("0" <= ch <= "9") or ch in "-_.") else "_" for ch in raw
    ).strip("._")
    return cleaned or "default"


def _windows_user_segment(environ: Optional[Mapping[str, str]] = None) -> str:
    env = os.environ if environ is None else environ
    for key in _WINDOWS_USER_ENV_KEYS:
        username = env.get(key, "").strip()
        if username:
            return _sanitize_user_segment(username)
    return "default"


def _user_segment() -> str:
    if os.name != "nt" and hasattr(os, "getuid"):
        return str(os.getuid())
    return _windows_user_segment()


def default_runtime_dir(root: Optional[Path] = None) -> Path:
    """返回当前用户的共享 per-user runtime 目录路径（Detached Mode，不自动创建）。"""
    base = Path(root) if root is not None else _runtime_root()
    return base / f"{RUNTIME_DIR_PREFIX}{_user_segment()}"


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


def acquire_runtime(root: Optional[Path] = None, *, shared: bool = True) -> RuntimeLease:
    """获取 runtime 目录并校验安全性。

    - ``shared=True``（Detached Mode）：获取或创建共享 per-user 目录（``swanlab-core-<uid>``），
      供同节点同用户的多 Rank / 多 Run 共享复用；
    - ``shared=False``（Owner Mode）：在 ``runtime_root`` 下创建每 Run 独占的 0700 临时目录，
      避免与共享 Core 或其他 Owner Run 冲突，退出时随 Run 清理。

    UDS 路径超长时直接失败，不回退 cwd。
    """
    base = Path(root) if root is not None else _runtime_root()
    created = False
    if shared:
        runtime_dir = default_runtime_dir(base)
        try:
            os.mkdir(runtime_dir, RUNTIME_DIR_MODE)
            created = True
        except FileExistsError:
            pass
        except OSError as exc:
            raise _runtime_error(f"cannot create runtime dir {runtime_dir}: {exc}", cause=exc) from exc
    else:
        try:
            runtime_dir = Path(tempfile.mkdtemp(prefix=f"{RUNTIME_DIR_PREFIX}own-", dir=str(base)))
            created = True
        except OSError as exc:
            raise _runtime_error(f"cannot create owner runtime dir under {base}: {exc}", cause=exc) from exc

    try:
        if created and os.name != "nt":
            # 消除 umask 影响，仅在本次新建时修正权限；已存在目录直接交由校验拒绝。
            os.chmod(runtime_dir, RUNTIME_DIR_MODE)
        _validate_runtime_dir(runtime_dir)
        if os.name != "nt":
            socket_path = runtime_dir / SOCKET_NAME
            if len(os.fsencode(str(socket_path))) > MAX_UNIX_PATH_BYTES:
                raise _runtime_error(f"unix socket path exceeds {MAX_UNIX_PATH_BYTES} bytes: {socket_path}")
    except BaseException:
        if created:
            shutil.rmtree(runtime_dir, ignore_errors=True)
        raise

    return RuntimeLease(
        dir=runtime_dir,
        port_file=runtime_dir / PORT_FILE_NAME,
        stderr_file=runtime_dir / STDERR_FILE_NAME,
        owner_pid=fork.current_pid(),
    )


def cleanup_runtime(lease: RuntimeLease) -> None:
    """清理隔离的 Owner / 临时 runtime 目录；共享 per-user 目录保留不删。"""
    try:
        if fork.is_forked(lease.owner_pid):
            return
        runtime_dir = Path(lease.dir)
        if not runtime_dir.name.startswith(RUNTIME_DIR_PREFIX):
            return
        if runtime_dir.name == f"{RUNTIME_DIR_PREFIX}{_user_segment()}":
            return
        if runtime_dir.is_symlink() or not runtime_dir.is_dir():
            return
        if hasattr(os, "getuid") and runtime_dir.lstat().st_uid != os.getuid():
            return
        shutil.rmtree(runtime_dir)
    except Exception as exc:  # noqa: BLE001 - 清理失败不能覆盖原始异常
        console.warning(f"failed to cleanup core runtime {lease.dir}: {exc}")


def try_acquire_start_lock(
    runtime_dir: Path,
    *,
    now_ns_fn: Optional[Callable[[], int]] = None,
) -> Optional[StartupLock]:
    """在 runtime 目录内尝试非阻塞获取启动锁；被占用时返回 None。"""
    return _try_acquire_start_lock(
        runtime_dir,
        validate_dir=_validate_runtime_dir,
        now_ns_fn=now_ns_fn,
    )


def acquire_start_lock(
    runtime_dir: Path,
    *,
    timeout: float = START_LOCK_TIMEOUT,
    poll_interval: float = START_LOCK_POLL_INTERVAL,
    now_ns_fn: Optional[Callable[[], int]] = None,
    sleep_fn: Optional[Callable[[float], None]] = None,
) -> StartupLock:
    """在 runtime 目录内竞争启动锁，直到成功或超时抛出 CoreGoError。"""
    return _acquire_start_lock(
        runtime_dir,
        validate_dir=_validate_runtime_dir,
        timeout=timeout,
        poll_interval=poll_interval,
        now_ns_fn=now_ns_fn,
        sleep_fn=sleep_fn,
    )
