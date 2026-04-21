from functools import wraps
from typing import Dict, List, Optional, Set, Tuple

from swanlab.api.typings.common import ApiIdentityEnum


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


def parse_column_type(column: str) -> str:
    """从前缀中获取指标类型"""
    column_type = column.split(".", 1)[0]
    if column_type == "summary":
        return "SCALAR"
    elif column_type == "config":
        return "CONFIG"
    else:
        return "STABLE"


def to_camel_case(name: str) -> str:
    """将下划线命名转化为驼峰命名"""
    return "".join([w.capitalize() if i > 0 else w for i, w in enumerate(name.split("_"))])


#TODO: 私有化接口装饰器
# def with_self_hosted(identity: ApiIdentityEnum = "user"):
#     """
#     用于需要在私有化环境下使用的接口的装饰器。
#     :param identity: 用户身份，默认为 "user"，如果为 "root"，则会额外验证是否为根用户。
#     """

#     def decorator(func):
#         @wraps(func)
#         def wrapper(self, *args, **kwargs):
#             client = getattr(self, "_client", None)
#             if not isinstance(client, Client):
#                 raise AttributeError("There is no SwanLab client instance.")

#             # 1. 尝试获取私有化服务信息
#             try:
#                 self_hosted_info = get_self_hosted_init(client)
#             except ApiError:
#                 raise ValueError("You haven't launched a swanlab self-hosted instance. This usages are not available.")

#             if not self_hosted_info.get("enabled", False):
#                 raise ValueError("SwanLab self-hosted instance hasn't been ready yet.")
#             if self_hosted_info.get("expired", True):
#                 raise ValueError("SwanLab self-hosted instance has expired.")

#             # 2. 检测用户权限（商业版root用户功能）
#             if identity == "root":
#                 if not self_hosted_info.get("root", False):
#                     raise ValueError("You don't have permission to perform this action. Please login as a root user")
#                 if not getattr(self, "is_self", True):
#                     raise ValueError("This root-only action can only be performed by the logged-in root user.")

#             return func(self, *args, **kwargs)

#         return wrapper

#     return decorator


_SPECIAL_FILTER_MAP = {
    # (backend_key, operator) — 用户侧 key 到后端字段名和操作符的映射
    # backend_key: 后端 API 实际接受的字段名
    # operator: 筛选操作符，EQ=精确匹配，IN=包含匹配（用于 tags 列表）
    "group": ("cluster", "EQ"),
    "tags": ("labels", "IN"),
    "name": ("name", "EQ"),
    "username": ("user.username", "EQ"),
    "job_type": ("job", "EQ"),
}


def parse_filter(key: str, value: object) -> Dict[str, object]:
    """将用户侧筛选条件转换为后端 filter 格式。

    :param key: 筛选字段名。预定义字段（group/tags/name/username/job_type）会映射到后端字段名；
        其他字段按 column type 自动转换：STABLE 类型转 camelCase，其余取最后一段。
    :param value: 筛选值。预定义字段中 tags 接受列表/元组，其余均为单值（内部统一包装为列表）。
    :return: 后端 filter 字典，包含 key / active / value / op / type 五个字段。
    """
    if key in _SPECIAL_FILTER_MAP:
        backend_key, op = _SPECIAL_FILTER_MAP[key]
        filter_value = list(value) if key == "tags" and isinstance(value, (list, tuple)) else [value]
        return {"key": backend_key, "active": True, "value": filter_value, "op": op, "type": "STABLE"}
    ct = parse_column_type(key)
    return {
        "key": to_camel_case(key) if ct == "STABLE" else key.split(".", 1)[-1],
        "active": True,
        "value": [value],
        "op": "EQ",
        "type": ct,
    }


def unwrap_api_payload(data):
    """提取 raw resp 的 data 响应."""
    if isinstance(data, dict) and "data" in data and isinstance(data["data"], (dict, list)):
        return data["data"]
    return data


# mulitpart-save
def extract_upload_id(payload: Dict[str, object]) -> Optional[str]:
    upload_id = payload.get("uploadId")
    if isinstance(upload_id, str) and upload_id:
        return upload_id
    return None


# multipart-save
def extract_part_urls(payload: Dict[str, object]) -> List[Tuple[int, str]]:
    parts = payload.get("parts")
    if not isinstance(parts, list):
        raise ValueError("Multipart upload URLs are missing in prepare response.")

    resolved = []
    for part in parts:
        if not isinstance(part, dict):
            raise ValueError("Multipart prepare response contains invalid part data.")
        part_number = part.get("partNumber")
        url = part.get("url")
        if not isinstance(part_number, int) or not isinstance(url, str) or not url:
            raise ValueError("Invalid partNumber or url in multipart response.")
        resolved.append((part_number, url))

    return sorted(resolved, key=lambda item: item[0])
