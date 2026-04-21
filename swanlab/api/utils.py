"""
@author: caddiesnew
@file: utils.py
@time: 2026/4/20
@description: 公共查询 API 工具函数
"""

from typing import Dict, NamedTuple, Optional, Set


class Label(NamedTuple):
    name: str

    def __str__(self) -> str:
        return self.name


def get_properties(obj: object, _visited: Optional[Set[int]] = None) -> Dict[str, object]:
    """递归获取实例中所有 property 的值，用于 to_dict() 默认实现。"""
    if _visited is None:
        _visited = set()
    obj_id = id(obj)
    if obj_id in _visited:
        return {}
    _visited = _visited | {obj_id}

    result = {}
    for name in dir(obj):
        if name.startswith("_"):
            continue
        if isinstance(getattr(type(obj), name, None), property):
            value = getattr(obj, name, None)
            result[name] = value if type(value).__module__ == "builtins" else get_properties(value, _visited)
    return result
