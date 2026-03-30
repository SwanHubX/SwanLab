"""
@author: cunyue
@file: conda.py
@time: 2026/3/30
@description: Conda 环境信息采集
"""

import subprocess
from typing import Optional


def get_metadata_conda() -> Optional[str]:
    """获取 conda 环境信息"""
    result = subprocess.run(["conda", "env", "export"], capture_output=True, text=True, timeout=15, check=True)
    return result.stdout
