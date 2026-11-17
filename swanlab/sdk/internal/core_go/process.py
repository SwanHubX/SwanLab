"""
@author: Nexisato
@file: process.py
@time: 2026/10/27 15:40
@description: swanlab-core OS 子进程生命周期管理。

负责 Popen 拉起、等待 port-file 就绪、早死分类、有界 stderr 诊断以及
wait → terminate → kill 的分级回收。不解析业务配置，不直接 import grpc。
"""

import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.core_go.config import (
    ENDPOINT_POLL_INTERVAL,
    EXIT_ENDPOINT_BUSY,
    EXIT_USAGE_ERROR,
    KILL_TIMEOUT,
    PORT_FILE_TIMEOUT,
    STDERR_TAIL_BYTES,
    TERMINATE_TIMEOUT,
    WAIT_TIMEOUT,
    CoreGoError,
    ReasonCode,
    Stage,
)
from swanlab.sdk.internal.core_go.runtime import (
    Endpoint,
    RuntimeLease,
    read_port_file,
    validate_core_binary,
    validate_ownership,
)
from swanlab.sdk.internal.pkg import fork

__all__ = ["reap_process", "spawn_core", "stderr_tail", "wait_endpoint"]


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


def _format_go_duration(seconds: float) -> str:
    """将非负秒数格式化为 Go time.ParseDuration 可解析的固定小数/毫秒串，避免科学计数法。"""
    if seconds <= 0:
        return "0s"
    millis = max(1, round(seconds * 1000))
    if millis % 1000 == 0:
        return f"{millis // 1000}s"
    return f"{millis}ms"


def spawn_core(
    binary: Path,
    lease: RuntimeLease,
    parent_pid: Optional[int] = None,
    *,
    detach: bool = False,
    idle_timeout: float = 0.0,
) -> subprocess.Popen:
    """拉起 swanlab-core，返回 Popen。

    Owner Mode 下只传 --port-filename / --parent-pid；Attach Mode（detach=True）
    传 --port-filename / --detach（及正数 --idle-timeout），不绑定父进程生命周期。
    子进程 cwd 是 runtime 目录，因此相对路径必须先按调用方 cwd 解析为绝对路径。
    """
    if idle_timeout < 0:
        raise ValueError("idle_timeout must be non-negative")
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
    ]
    if detach:
        argv.append("--detach")
        if idle_timeout > 0:
            argv.extend(["--idle-timeout", _format_go_duration(idle_timeout)])
    else:
        argv.extend(["--parent-pid", str(parent_pid if parent_pid is not None else fork.current_pid())])

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


def _early_exit_error(lease: RuntimeLease, exit_code: Optional[int]) -> CoreGoError:
    if exit_code == EXIT_ENDPOINT_BUSY:
        reason = ReasonCode.ENDPOINT_BUSY
        message = "core exited because the endpoint is already in use"
    elif exit_code == EXIT_USAGE_ERROR:
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
            raise _early_exit_error(lease, exit_code)

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
                raise _early_exit_error(lease, exit_code)
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
