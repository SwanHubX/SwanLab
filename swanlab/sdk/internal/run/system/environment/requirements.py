"""
@author: cunyue
@file: requirements.py
@time: 2026/3/30
@description: 依赖包信息采集
"""

import subprocess
from typing import Optional

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.utils.helper import catch_and_return_none


def _try_run(cmd: list, timeout: int = 5, check: bool = False) -> Optional[subprocess.CompletedProcess]:
    """尝试执行命令，命令不存在时返回 None 而不是抛异常"""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=check)
    except FileNotFoundError:
        return None


@catch_and_return_none(on_error=lambda e: console.debug(f"Failed to get environment requirements: {e}"))
def get() -> str:
    """获取当前环境依赖"""
    # 尝试 pixi
    result = _try_run(["pixi", "list"])
    if result and result.returncode == 0:
        return result.stdout

    # 尝试 uv
    result = _try_run(["uv", "pip", "list", "--format=freeze"])
    if result and result.returncode == 0:
        return result.stdout

    # 尝试 pip
    result = _try_run(["pip", "list", "--format=freeze"], timeout=15, check=True)
    if result:
        return result.stdout

    raise ValueError("Failed to get environment requirements")
