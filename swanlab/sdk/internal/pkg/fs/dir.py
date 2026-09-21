"""
@author: cunyue
@file: dir.py
@time: 2026/3/11 13:44
@description: SwanLab SDK 目录辅助函数
"""

import os
import tempfile
import time
from pathlib import Path
from typing import Union

from swanlab.sdk.internal.pkg import console, safe


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


def _probe_writable(p: Path) -> None:
    """
    目录可写性探针：创建命名文件 → 写入 → 关闭 → 删除。

    不使用 tempfile.TemporaryFile：其匿名文件（O_TMPFILE）或「fd 仍打开时立即
    unlink」的语义在部分 NAS 上会返回 EIO。只有创建、写入和关闭
    用于判定目录可写；删除仅尽力清理本次探针文件。

    实现思路：
    1. 用 mkstemp 在目标目录创建命名探针文件。
    2. 写入一个字节，验证文件可写。
    3. 先关闭 fd，避免部分 FUSE / NAS 在文件仍打开时 unlink 返回 EIO。
    4. 只尝试删除本次探针文件，不扫描或清理历史文件。
    5. 创建、写入或关闭失败时上抛供 safe_mkdir 重试；删除失败仅记录 trace。
    """
    fd, name = tempfile.mkstemp(dir=p, prefix=PROBE_PREFIX)
    try:
        try:
            os.write(fd, b"0")
        except BaseException:
            # close 的异常不得覆盖原始错误
            with safe.block(OSError, message=None):
                os.close(fd)
            raise
        os.close(fd)
    finally:
        # 目录已经通过创建、写入和关闭完成可写性探测。
        # 清理失败不参与可写性判定，也不扫描或删除历史探针文件。
        with safe.block(
            OSError,
            message=f"Failed to clean up writability probe file [{name}]",
            write_to_tty=False,
        ):
            os.unlink(name)


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
    while True:
        try:
            _probe_writable(p)
            break
        except PermissionError as e:
            # 创建、写入或关闭探针文件被拒：权限不足不会因重试而恢复。
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
