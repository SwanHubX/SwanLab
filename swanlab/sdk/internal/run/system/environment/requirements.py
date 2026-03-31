"""
@author: cunyue
@file: requirements.py
@time: 2026/3/30
@description: 依赖包信息采集
"""

import subprocess

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.utils.helper import catch_and_return_none


@catch_and_return_none(on_error=lambda e: console.debug(f"Failed to get environment requirements: {e}"))
def get() -> str:
    """获取当前环境依赖"""
    # 尝试 pixi
    result = subprocess.run(["pixi", "list"], capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        return result.stdout

    # 尝试 uv
    result = subprocess.run(["uv", "pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        return result.stdout

    # 尝试 pip
    result = subprocess.run(["pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=15, check=True)
    return result.stdout
