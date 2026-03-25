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
    column_type = column.split('.', 1)[0]
    if column_type == 'summary':
        return 'SCALAR'
    elif column_type == 'config':
        return 'CONFIG'
    else:
        return 'STABLE'


# 将下划线命名转化为驼峰命名
def to_camel_case(name: str) -> str:
    return ''.join([w.capitalize() if i > 0 else w for i, w in enumerate(name.split('_'))])


def unwrap_api_payload(data):
    if isinstance(data, dict) and "data" in data and isinstance(data["data"], (dict, list)):
        return data["data"]
    return data


def extract_file_payloads(payload) -> List[Dict[str, object]]:
    if isinstance(payload, list):
        if not all(isinstance(item, dict) for item in payload):
            raise ValueError("Unexpected file payload shape from save API.")
        return payload
    if isinstance(payload, dict):
        if any(isinstance(payload.get(key), str) for key in ("uploadUrl", "upload_url", "url", "presignedUrl")):
            return [payload]
        if any(key in payload for key in ("uploadId", "upload_id", "parts", "uploadUrls", "upload_urls", "urls")):
            return [payload]
        for key in ("files", "items", "list"):
            value = payload.get(key)
            if isinstance(value, list) and all(isinstance(item, dict) for item in value):
                return value
    raise ValueError("Unexpected file payload shape from save API.")


def extract_upload_url(payload: Dict[str, object]) -> str:
    for key in ("uploadUrl", "upload_url", "url", "presignedUrl", "presigned_url"):
        value = payload.get(key)
        if isinstance(value, str) and value != "":
            return value
    raise ValueError("Upload URL is missing in prepare response.")


def extract_upload_id(payload: Dict[str, object]) -> Optional[str]:
    for key in ("uploadId", "upload_id"):
        value = payload.get(key)
        if isinstance(value, str) and value != "":
            return value
    return None


def extract_part_size(payload: Dict[str, object], default_size: int) -> int:
    for key in ("partSize", "part_size"):
        value = payload.get(key)
        if isinstance(value, int) and value > 0:
            return value
    return default_size


def extract_part_urls(payload: Dict[str, object]) -> List[Tuple[int, str]]:
    parts = payload.get("parts")
    if isinstance(parts, list):
        resolved = []
        for index, part in enumerate(parts, start=1):
            if not isinstance(part, dict):
                raise ValueError("Multipart prepare response contains invalid part data.")
            number = part.get("partNumber", part.get("part_number", index))
            resolved.append((int(number), extract_upload_url(part)))
        return sorted(resolved, key=lambda item: item[0])

    urls = payload.get("uploadUrls", payload.get("upload_urls", payload.get("urls")))
    if isinstance(urls, list):
        resolved = []
        for index, url in enumerate(urls, start=1):
            if not isinstance(url, str) or url == "":
                raise ValueError("Multipart prepare response contains invalid upload URL.")
            resolved.append((index, url))
        return resolved

    raise ValueError("Multipart upload URLs are missing in prepare response.")


__all__ = [
    "parse_column_type",
    "to_camel_case",
    "unwrap_api_payload",
    "extract_file_payloads",
    "extract_upload_url",
    "extract_upload_id",
    "extract_part_size",
    "extract_part_urls",
]
