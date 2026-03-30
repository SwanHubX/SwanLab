"""
@author: cunyue
@file: conda.py
@time: 2026/3/30
@description: Conda 环境信息采集
"""

import subprocess

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.utils.helper import catch_and_return_none


@catch_and_return_none(on_error=lambda e: console.error("Failed to get conda environment information: %s", str(e)))
def get() -> str:
    """获取 conda 环境信息"""
    result = subprocess.run(["conda", "env", "export"], capture_output=True, text=True, timeout=15, check=True)
    return result.stdout
