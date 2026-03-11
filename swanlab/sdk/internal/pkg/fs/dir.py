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


def safe_mkdirs(*paths: Union[str, Path], timeout: float = TIMEOUT):
    """
    安全地创建多个目录。
    :param paths: 目录路径列表
    :param timeout: 超时时间（秒）
    """
    for path in paths:
        safe_mkdir(path, timeout=timeout)


def safe_mkdir(path: Union[str, Path], timeout: float = TIMEOUT) -> Path:
    """
    安全地创建目录，带有抗 NAS 异步延迟的探针机制。
    :param path: 目录路径
    :param timeout: 超时时间（秒）
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    start_time = time.time()

    while not (p.exists() and p.is_dir()):
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Directory creation timeout (NAS delay too high): {p}")
        time.sleep(0.05)

    while True:
        try:
            with tempfile.TemporaryFile(dir=p, prefix=".swanlab_test_") as f:
                f.write(b"0")
            break
        except OSError as e:
            if time.time() - start_time > timeout:
                console.error(f"Directory {p} exists but is not writable yet. NAS sync issue? Error: {e}")
                raise TimeoutError(f"Directory {p} is not writable within {timeout}s.")
            time.sleep(0.1)

    return p
