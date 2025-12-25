"""
@author: Zhou Qiyang
@file: utils.py
@time: 2025/12/25 21:32
@description: OpenApi 中使用的工具函数
"""

from typing import Dict, List


def flatten_runs(runs: Dict) -> List:
    """
    展开分组后的实验数据，返回一个包含所有实验的列表
    """
    flat_runs = []
    for group in runs.values():
        if isinstance(group, Dict):
            flat_runs.extend(flatten_runs(group))
        else:
            flat_runs.extend(group)
    return flat_runs
