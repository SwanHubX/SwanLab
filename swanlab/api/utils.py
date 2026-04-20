"""
@author: caddiesnew
@file: utils.py
@time: 2026/4/20
@description: 公共查询 API 工具函数
"""

from dataclasses import dataclass
from typing import Dict


@dataclass
class Label:
    name: str

    def __str__(self) -> str:
        return self.name


def get_properties(obj: object) -> Dict[str, object]:
    """递归获取实例中所有 property 的值，用于 to_dict() 默认实现。"""
    result = {}
    for name in dir(obj):
        if name.startswith("_"):
            continue
        if isinstance(getattr(type(obj), name, None), property):
            value = getattr(obj, name, None)
            result[name] = value if type(value).__module__ == "builtins" else get_properties(value)
    return result


def parse_column_type(column: str) -> str:
    """从前缀中获取指标列类型。"""
    column_type = column.split(".", 1)[0]
    if column_type == "summary":
        return "SCALAR"
    elif column_type == "config":
        return "CONFIG"
    else:
        return "STABLE"


def to_camel_case(name: str) -> str:
    """将下划线命名转化为驼峰命名。"""
    return "".join(w.capitalize() if i > 0 else w for i, w in enumerate(name.split("_")))
