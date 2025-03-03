"""
@author:    cunyue
@file:      conda.py
@time:      2025/3/2 00:07
@description: 收集conda信息
"""

import subprocess
from typing import Optional

import yaml


def get_conda() -> Optional[str]:
    """
    获取conda信息，如果不存在则返回None
    Returns: str
    """
    result = subprocess.run(["conda env export"], shell=True, capture_output=True, text=True, timeout=15)
    if result.returncode != 0:
        return None
    output = result.stdout
    if yaml.safe_load(output):
        return output
    return None
