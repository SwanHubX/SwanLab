"""
@author: cunyue
@file: conda.py
@time: 2026/3/30
@description: Conda 环境信息采集
"""

import subprocess

from swanlab.sdk.internal.pkg.safe import safe


@safe(level="debug", message="Failed to get conda environment")
def get() -> str:
    """获取 conda 环境信息"""
    result = subprocess.run(["conda", "env", "export"], capture_output=True, text=True, timeout=15, check=True)
    return result.stdout
