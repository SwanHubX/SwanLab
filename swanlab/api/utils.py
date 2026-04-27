"""
@author: caddiesnew
@file: utils.py
@time: 2026/4/20
@description: swanlab/api 实体层工具函数
"""

import re
from typing import Any, Dict, List, Optional, Set, Tuple, Type, get_args, get_type_hints

from swanlab.api.typings.common import (
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiColumnScalarTypeLiteral,
    ApiFilterOpLiteral,
    ApiFilterStableKeyLiteral,
    ApiMetricAllTypeLiteral,
    ApiMetricLogLevelLiteral,
    ApiSidebarLiteral,
    ApiSortOrderLiteral,
    ApiVisibilityLiteral,
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


# 路径解析
def resolve_run_path(path: str) -> Tuple[str, str]:
    """ "path like: user/proj_name/run_slug"""
    proj_path, run_slug = "", ""
    parts = path.split("/")
    if len(parts) != 3:
        return proj_path, run_slug
    run_slug = parts[-1]
    proj_path = path.rsplit("/", 1)[0]
    return (
        proj_path,
        run_slug,
    )


def validate_api_path(path: str, *, segments: int, label: str) -> None:
    """校验公开 API 入口 path 参数的段数和空白字符。"""
    if not isinstance(path, str):
        raise ValueError(f"{label} path must be a string")
    parts = path.split("/")
    if path != path.strip() or len(parts) != segments or any(part != part.strip() or not part for part in parts):
        raise ValueError(f"{label} path must contain {segments} non-empty segment(s), got {path!r}")


def validate_non_empty_string(value: str, *, label: str) -> None:
    """校验公开 API 入口中的非空字符串参数。"""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")


# ---------------------------------------------------------------------------
# POST /runs/shows 参数校验常量（从 typings 中的 Literal 类型提取，避免重复定义）
# ---------------------------------------------------------------------------

_VALID_SIDEBAR_TYPES = frozenset(get_args(ApiSidebarLiteral))
_VALID_OPS = frozenset(get_args(ApiFilterOpLiteral))
_VALID_ORDERS = frozenset(get_args(ApiSortOrderLiteral))
_STABLE_KEYS = frozenset(get_args(ApiFilterStableKeyLiteral))
_VALID_VISIBILITIES = frozenset(get_args(ApiVisibilityLiteral))

_PROJECT_NAME_RE = re.compile(r"^[0-9a-zA-Z\-_.+]+$")

# 列相关校验常量
_VALID_COLUMN_CLASSES = frozenset(get_args(ApiColumnClassLiteral))
_VALID_COLUMN_DATA_TYPES = frozenset(get_args(ApiColumnDataTypeLiteral))
_VALID_COLUMN_SCALAR_TYPES = frozenset(get_args(ApiColumnScalarTypeLiteral))

# 指标相关校验常量
_VALID_METRIC_ALL_TYPES = frozenset(get_args(ApiMetricAllTypeLiteral))
_VALID_METRIC_LOG_LEVELS = frozenset(get_args(ApiMetricLogLevelLiteral))


def _check_required(item: Dict[str, Any], keys: Set[str]) -> None:
    if not isinstance(item, dict):
        raise ValueError(f"Expected dict item, got {type(item).__name__}")
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


def validate_metric_type(metric_type: str, key: Optional[str] = None) -> None:
    """校验 metric_type 的合法性。非 LOG 类型必须提供非空 key。"""
    if metric_type not in _VALID_METRIC_ALL_TYPES:
        raise ValueError(f"Invalid metric_type: {metric_type!r}, expected one of {sorted(_VALID_METRIC_ALL_TYPES)}")
    if metric_type != "LOG" and (not isinstance(key, str) or not key.strip()):
        raise ValueError(f"key is required for metric_type {metric_type!r}, got key={key!r}")


def validate_metric_log_level(level: str) -> None:
    """校验 metric log level 的合法性。"""
    if level not in _VALID_METRIC_LOG_LEVELS:
        raise ValueError(f"Invalid metric log level: {level!r}, expected one of {sorted(_VALID_METRIC_LOG_LEVELS)}")


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


def validate_update_active(
    items: Optional[List[Dict[str, Any]]],
    validator,
    *,
    label: str = "items",
) -> List[Dict[str, Any]]:
    """校验每个 item 并补充 active: True，返回可直接发送的列表。"""
    if items is None:
        return []
    if not isinstance(items, list):
        raise ValueError(f"{label} must be a list")
    if not items:
        return []
    for item in items:
        if not isinstance(item, dict):
            raise ValueError(f"{label} items must be dicts")
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


def parse_column_data_type(column_type: str):
    """解析列类型。"""
    validate_column_params(column_type=column_type)
    if column_type in _VALID_COLUMN_SCALAR_TYPES:
        return "SCALAR"
    # 新加入的类型默认指定为 media
    return "MEDIA"


# ---------------------------------------------------------------------------
# 创建项目 / 实验的参数校验
# ---------------------------------------------------------------------------


def validate_project_name(name: str) -> None:
    if not 1 <= len(name) <= 100:
        raise ValueError("Project name must be between 1 and 100 characters.")
    if not _PROJECT_NAME_RE.match(name):
        raise ValueError("Project name can only contain 0-9, a-z, A-Z, -, _, ., +")


def validate_visibility(visibility: str) -> None:
    """校验 visibility 的合法性。"""
    if visibility not in _VALID_VISIBILITIES:
        raise ValueError(f"Invalid visibility: {visibility!r}, expected one of {sorted(_VALID_VISIBILITIES)}")


def validate_metric_keys(keys: List[str]) -> None:
    """校验 metric keys 列表的合法性。"""
    if not isinstance(keys, list) or not keys or any(not isinstance(key, str) or not key.strip() for key in keys):
        raise ValueError("keys must be a non-empty list of non-empty strings")
