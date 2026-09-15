"""
@author: cunyue
@file: dir.py
@time: 2026/3/11 13:44
@description: SwanLab SDK 目录辅助函数
"""

import contextlib
import os
import tempfile
import time
from pathlib import Path
from typing import List, Union

from .. import console


def _get_fs_timeout(default: float = 5.0) -> float:
    """
    安全地从环境变量获取文件系统超时时间
    防范用户输入非数字 ("abc") 或非法数字 (-1.0) 导致模块导入崩溃
    """
    env_val = os.environ.get("SWANLAB_FS_TIMEOUT")
    if env_val is not None:
        try:
            val = float(env_val)
            if val > 0:
                return val
            else:
                # 用户如果设了负数或0，打个警告，回退到默认值
                console.warning(f"SWANLAB_FS_TIMEOUT must be > 0, got {val}. Using default {default}s.")
        except ValueError:
            console.warning(f"Invalid SWANLAB_FS_TIMEOUT value: '{env_val}'. Using default {default}s.")
    return default


# 模块加载时安全获取
TIMEOUT = _get_fs_timeout()

# 可写性探针临时文件前缀
PROBE_PREFIX = ".swanlab_test_"


def safe_mkdirs(*paths: Union[str, Path], timeout: float = TIMEOUT, ensure_clean: bool = False):
    """
    安全地创建多个目录。
    :param paths: 目录路径列表
    :param timeout: 超时时间（秒）
    :param ensure_clean: 如果为 True，要求目标目录在创建前不存在
    """
    for path in paths:
        safe_mkdir(path, timeout=timeout, ensure_clean=ensure_clean)


def _probe_writable(p: Path, stale: List[str]) -> None:
    """
    目录可写性探针：创建命名文件 → 写入 → 关闭 → 删除。

    不使用 tempfile.TemporaryFile：其匿名文件（O_TMPFILE）或「fd 仍打开时立即
    unlink」的语义在部分 NAS 上会返回 EIO。清理动作不掩盖原始异常。
    """
    if stale:
        # 上一轮自己的残留：先删（FileNotFoundError = 已被清掉，视为干净）
        try:
            os.unlink(stale[0])
        except FileNotFoundError:
            pass
        stale.clear()
        # 其他 OSError 如实抛出，由上层重试

    fd, name = tempfile.mkstemp(dir=p, prefix=PROBE_PREFIX)
    try:
        try:
            os.write(fd, b"0")
        except BaseException:
            # close 的异常不得覆盖原始错误
            with contextlib.suppress(OSError):
                os.close(fd)
            raise
        os.close(fd)
    except BaseException:
        # 任意失败（含 KeyboardInterrupt）：尽力清理，不掩盖原始异常
        try:
            os.unlink(name)
        except OSError:
            stale.append(name)
        raise
    # 删除失败如实上报（文件名留在 stale 里），由上层重试
    try:
        os.unlink(name)
    except OSError:
        stale.append(name)
        raise


def safe_mkdir(path: Union[str, Path], timeout: float = TIMEOUT, ensure_clean: bool = False) -> Path:
    """
    安全地创建目录，带有抗异步文件系统延迟的探针机制。

    创建后会探测目录是否真正可见且可写，以容忍 NAS / NFS 等的异步IO延迟，可写性探测使用命名临时文件。
    权限不足属于不可恢复错误，会立即抛出 PermissionError，不会重试。

    :param path: 目录路径
    :param timeout: 超时时间（秒）
    :param ensure_clean: 如果为 True，要求目标目录在创建前不存在（原子检查，通过 mkdir(exist_ok=False) 实现）
    :raises FileExistsError: ensure_clean=True 且目录已存在时抛出
    :raises PermissionError: 目录不可创建或不可写时抛出
    """
    p = Path(path)
    try:
        if ensure_clean:
            p.mkdir(parents=True, exist_ok=False)
        else:
            p.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        raise PermissionError(
            f"Cannot create directory [{p}]: permission denied. "
            "Please choose a writable log_dir or update directory permissions."
        ) from None

    start_time = time.time()

    # 探测一：目录创建后可能因异步延迟暂不可见
    while not (p.exists() and p.is_dir()):
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Directory creation timed out, filesystem may be slow or remote: {p}")
        time.sleep(0.05)

    # 探测二：目录可见后可能仍暂不可写（权限问题立即失败，其余重试到超时）
    stale: List[str] = []  # 上次尝试未删掉的探针文件名（如有），重试时只精确清理这一个
    while True:
        try:
            _probe_writable(p, stale)
            break
        except PermissionError as e:
            # 只可能是探针本体（mkstemp/写/删自己的文件）被拒：权限不足不会因
            # 重试而恢复，立即失败。探针从不 unlink 其他进程文件，无清理阶段误伤。
            console.trace(f"Directory [{p}] is not writable, underlying error: {e}")
            raise PermissionError(
                f"Directory [{p}] is not writable. Please choose a writable log_dir or update directory permissions."
            ) from None
        except OSError as last_error:
            # 取最后一次原始 OSError，避免把硬错误误报成单纯超时
            if time.time() - start_time > timeout:
                console.trace(f"Directory {p} exists but is not writable within {timeout}s")
                raise TimeoutError(
                    f"Directory [{p}] exists but is not writable within {timeout}s, FILESYSTEM may be slow or remote."
                ) from last_error
            time.sleep(0.1)

    return p
