"""
@author: caddiesnew
@file: utils.py
@time: 2026/4/20
@description: swanlab/api 实体层工具函数
"""

from typing import Any, Dict, List, Optional, Set, Type, get_args, get_type_hints

from swanlab.api.typings.common import ApiColumnClassLiteral, ApiColumnDataTypeLiteral, ApiSidebarLiteral
from swanlab.api.typings.experiment import (
    ApiFilterOpLiteral,
    ApiSortOrderLiteral,
    ApiStableKeyLiteral,
)


def strip_dict(data: Any, typed_cls: Type) -> Dict[str, Any]:
    """将原始 API 响应字典裁剪为只保留 TypedDict 中声明的字段。"""
    if not data:
        return {}
    hints = get_type_hints(typed_cls) if hasattr(typed_cls, "__annotations__") else {}
    return {k: data[k] for k in hints if k in data}


def get_properties(obj: object, _visited: Optional[Set[int]] = None) -> Dict[str, object]:
    """递归获取实例中所有 property 的值，用于 json() 默认实现。"""
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


# ---------------------------------------------------------------------------
# POST /runs/shows 参数校验常量（从 typings 中的 Literal 类型提取，避免重复定义）
# ---------------------------------------------------------------------------

_VALID_SIDEBAR_TYPES = frozenset(get_args(ApiSidebarLiteral))
_VALID_OPS = frozenset(get_args(ApiFilterOpLiteral))
_VALID_ORDERS = frozenset(get_args(ApiSortOrderLiteral))
_STABLE_KEYS = frozenset(get_args(ApiStableKeyLiteral))

# 列相关校验常量
_VALID_COLUMN_CLASSES = frozenset(get_args(ApiColumnClassLiteral))
_VALID_COLUMN_DATA_TYPES = frozenset(get_args(ApiColumnDataTypeLiteral))


def _check_required(item: Dict[str, Any], keys: Set[str]) -> None:
    missing = keys - item.keys()
    if missing:
        raise ValueError(f"Missing required fields: {sorted(missing)}, got {sorted(item.keys())}")


def _check_type_field(item: Dict[str, Any]) -> None:
    t = item.get("type", "")
    if t not in _VALID_SIDEBAR_TYPES:
        raise ValueError(f"Invalid type: {t!r}, expected one of {sorted(_VALID_SIDEBAR_TYPES)}")


def _check_stable_key(item: Dict[str, Any]) -> None:
    if item.get("type") == "STABLE" and item["key"] not in _STABLE_KEYS:
        raise ValueError(f"Invalid STABLE key: {item['key']!r}, expected one of {sorted(_STABLE_KEYS)}")


def validate_filter(item: Dict[str, Any]) -> None:
    """校验单个 filter item 的合法性。"""
    _check_required(item, {"key", "type", "op", "value"})
    _check_type_field(item)
    _check_stable_key(item)
    if item["op"] not in _VALID_OPS:
        raise ValueError(f"Invalid filter op: {item['op']!r}, expected one of {sorted(_VALID_OPS)}")
    if not isinstance(item["value"], list):
        raise ValueError(f"filter value must be a list, got {type(item['value']).__name__}")


def validate_group(item: Dict[str, Any]) -> None:
    """校验单个 group item 的合法性。"""
    _check_required(item, {"key", "type"})
    _check_type_field(item)
    _check_stable_key(item)


def validate_sort(item: Dict[str, Any]) -> None:
    """校验单个 sort item 的合法性。"""
    _check_required(item, {"key", "type", "order"})
    _check_type_field(item)
    _check_stable_key(item)
    if item["order"] not in _VALID_ORDERS:
        raise ValueError(f"Invalid sort order: {item['order']!r}, expected one of {sorted(_VALID_ORDERS)}")


def _validate_and_build(
    items: Optional[List[Dict[str, Any]]],
    validator,
) -> List[Dict[str, Any]]:
    """校验每个 item 并补充 active: True，返回可直接发送的列表。"""
    if not items:
        return []
    for item in items:
        validator(item)
    return [{**item, "active": True} for item in items]


def validate_column_params(column_type: Optional[str] = None, column_class: Optional[str] = None) -> None:
    """
    校验列查询参数的合法性。

    :param column_type: 列的数据类型
    :param column_class: 列的分类
    :raises ValueError: 当参数不在允许的枚举值中时
    """
    if column_type is not None and column_type not in _VALID_COLUMN_DATA_TYPES:
        raise ValueError(f"Invalid column_type: {column_type!r}, expected one of {sorted(_VALID_COLUMN_DATA_TYPES)}")
    if column_class is not None and column_class not in _VALID_COLUMN_CLASSES:
        raise ValueError(f"Invalid column_class: {column_class!r}, expected one of {sorted(_VALID_COLUMN_CLASSES)}")
