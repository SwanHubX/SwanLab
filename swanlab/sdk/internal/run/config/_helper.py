"""
@author: cunyue
@file: _helper.py
@time: 2026/3/14 16:12
@description: SwanLabConfig 辅助函数
"""


def revert_config(config: dict) -> dict:
    """将 {key: {value, desc, sort}} 格式还原为 {key: value}，按 sort 升序排列。"""
    items = [(k, v["value"], v.get("sort", 0)) for k, v in config.items() if isinstance(v, dict) and "value" in v]
    return {k: v for k, v, _ in sorted(items, key=lambda x: x[2])}
