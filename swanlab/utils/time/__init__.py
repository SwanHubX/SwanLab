"""
@author: cunyue
@file: __init__.py
@time: 2026/8/18
@description: 时间戳解析辅助函数，供 API 查询 lane 与运行时 SDK lane 共用
"""

from datetime import datetime, timezone
from typing import Union

__all__ = ["TIMESTAMP_S_MAX", "TIMESTAMP_S_MIN", "parse_timestamp_s"]

# 秒级 Unix 时间戳的合法区间（10 位）
TIMESTAMP_S_MIN = 1_000_000_000
TIMESTAMP_S_MAX = 9_999_999_999


def parse_timestamp_s(value: Union[int, str]) -> int:
    """将时间统一转换为秒级 Unix 时间戳（10 位），用于 House 查询接口的 ``createdAt`` 参数。

    支持格式：
    - int / 数字字符串：秒（10 位）或毫秒（13 位）级时间戳，自动归一化为秒
    - ISO 8601 字符串：如 ``"2024-08-01T00:00:00Z"``、``"2024-08-01T08:00:00+08:00"``，无时区时按 UTC 解析

    :raises ValueError: ``value`` 为 None、空字符串、无法解析、非正数或越出合法区间时抛出；
        不做默认值回退，避免 House 退化为全历史扫描（慢查询）
    """
    if value is None:
        raise ValueError("timestamp value is required, got None")
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError("timestamp value must not be an empty string")
        if not text.isdigit():
            try:
                dt = datetime.fromisoformat(text.replace("Z", "+00:00"))
            except ValueError as exc:
                raise ValueError(f"Invalid ISO 8601 timestamp: {text!r}") from exc
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            value = int(dt.timestamp())
        else:
            value = int(text)
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"Expected positive int timestamp, got {value!r}")
    while value > TIMESTAMP_S_MAX:
        value //= 1000
    if not TIMESTAMP_S_MIN <= value <= TIMESTAMP_S_MAX:
        raise ValueError(f"Timestamp {value} out of range [{TIMESTAMP_S_MIN}, {TIMESTAMP_S_MAX}]")
    return value
