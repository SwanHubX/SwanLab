"""
@author: Zhou QiYang
@file: utils.py
@time: 2026/1/10 22:09
@description: 实验相关的后端API接口中的工具函数
"""

from typing import Dict, List, Optional, Tuple

from swanlab.core_python.api.type import ColumnType


# 从前缀中获取指标类型
def parse_column_type(column: str) -> ColumnType:
    column_type = column.split(".", 1)[0]
    if column_type == "summary":
        return "SCALAR"
    elif column_type == "config":
        return "CONFIG"
    else:
        return "STABLE"


# 将下划线命名转化为驼峰命名
def to_camel_case(name: str) -> str:
    return "".join(
        [w.capitalize() if i > 0 else w for i, w in enumerate(name.split("_"))]
    )


def unwrap_api_payload(data):
    if (
        isinstance(data, dict)
        and "data" in data
        and isinstance(data["data"], (dict, list))
    ):
        return data["data"]
    return data


def extract_upload_id(payload: Dict[str, object]) -> Optional[str]:
    upload_id = payload.get("uploadId")
    if isinstance(upload_id, str) and upload_id:
        return upload_id
    return None



def extract_part_urls(payload: Dict[str, object]) -> List[Tuple[int, str]]:
    urls = payload.get("urls")
    if not isinstance(urls, list):
        raise ValueError("Multipart upload URLs are missing in prepare response.")

    resolved = []
    for part in urls:
        if not isinstance(part, dict):
            raise ValueError("Multipart prepare response contains invalid part data.")
        part_number = part.get("partNumber")
        url = part.get("url")
        if not isinstance(part_number, int) or not isinstance(url, str) or not url:
            raise ValueError("Invalid partNumber or url in multipart response.")
        resolved.append((part_number, url))

    return sorted(resolved, key=lambda item: item[0])


__all__ = [
    "parse_column_type",
    "to_camel_case",
    "unwrap_api_payload",
    "extract_upload_id",
    "extract_part_urls",
]
