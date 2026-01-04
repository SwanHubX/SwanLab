"""
@author: Zhou QiYang
@file: utils.py
@time: 2026/1/4 18:03
@description: OpenApi 使用的常量和工具函数
"""

from typing import Dict, List

STATUS_OK = "OK"
STATUS_CREATED = "Created"


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
