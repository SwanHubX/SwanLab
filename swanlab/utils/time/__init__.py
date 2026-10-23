"""
@author: cunyue
@file: __init__.py
@time: 2026/8/18
@description: 时间戳解析辅助函数，供 API 查询 lane 与运行时 SDK lane 共用
"""

import math
import re
from datetime import datetime, timedelta, timezone, tzinfo
from typing import Optional, Union

__all__ = [
    "TIMESTAMP_S_MAX",
    "TIMESTAMP_S_MIN",
    "NANOSECONDS_PER_SECOND",
    "parse_timestamp_s",
    "parse_timestamp_ns",
    "parse_timezone",
]

# 秒级 Unix 时间戳的合法区间（10 位）
TIMESTAMP_S_MIN = 1_000_000_000
TIMESTAMP_S_MAX = 9_999_999_999

NANOSECONDS_PER_SECOND = 1_000_000_000

_UTC_NAMES = frozenset({"UTC", "Z", "GMT", "UCT", "UNIVERSAL"})
_OFFSET_RE = re.compile(r"^([+-])(\d{1,2})(?::?(\d{2}))?$")


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


def parse_timezone(value: Optional[Union[str, tzinfo]]) -> Optional[tzinfo]:
    """解析时区参数，返回 ``tzinfo``；``None`` 表示使用系统本地时区。

    支持格式：
    - ``None``：跟随系统本地时区
    - ``tzinfo`` 实例：原样返回
    - ``"UTC"`` / ``"Z"`` / ``"GMT"``：UTC
    - UTC 偏移：``"+08:00"``、``"-0530"``、``"+8"``
    - IANA 时区名：``"Asia/Shanghai"``（依赖系统时区数据库，Windows 无 tzdata 时请改用偏移）

    :raises ValueError: 时区为空、格式非法或 IANA 时区不存在时抛出
    """
    if value is None:
        return None
    if isinstance(value, tzinfo):
        return value
    if not isinstance(value, str):
        raise TypeError(f"Unsupported timezone type: {type(value).__name__}")
    text = value.strip()
    if not text:
        raise ValueError("timezone must not be an empty string")
    if text.upper() in _UTC_NAMES:
        return timezone.utc
    matched = _OFFSET_RE.match(text)
    if matched:
        sign = -1 if matched.group(1) == "-" else 1
        hours, minutes = int(matched.group(2)), int(matched.group(3) or 0)
        if hours > 23 or minutes > 59:
            raise ValueError(f"Invalid UTC offset: {text!r}")
        return timezone(timedelta(hours=hours, minutes=minutes) * sign)
    try:
        from zoneinfo import ZoneInfo, ZoneInfoNotFoundError
    except ImportError as exc:  # pragma: no cover - Python 3.9+ 自带 zoneinfo
        raise ValueError(
            f"Timezone database is unavailable, use a UTC offset like '+08:00' instead of {text!r}"
        ) from exc
    try:
        return ZoneInfo(text)
    except (ZoneInfoNotFoundError, ValueError, OSError) as exc:
        raise ValueError(
            f"Unknown timezone {text!r}. Use an IANA name like 'Asia/Shanghai' or a UTC offset like '+08:00'."
        ) from exc


def parse_timestamp_ns(value: Union[int, float, str, datetime], *, tz: Optional[tzinfo] = None) -> int:
    """将时间统一转换为纳秒级 Unix 时间戳，用于 ``swanlab sync`` 的时间范围过滤。

    支持格式：
    - int / float：秒级或毫秒级 Unix 时间戳，大于 10 位时自动按毫秒归一化为秒
    - 数字字符串：同上
    - ISO 8601 字符串：如 ``"2024-08-01T00:00:00Z"``、``"2024-08-01 08:00:00"``，也支持纯日期
    - ``datetime``：直接取时间戳

    不带时区的字符串和 ``datetime`` 按 *tz* 解析；*tz* 为 ``None`` 时跟随**系统本地时区**，
    与用户在终端输入墙上时间的习惯一致。带时区的输入始终以输入自带的时区为准。

    :param value: 待解析的时间
    :param tz: 无时区输入所使用的时区，``None`` 表示系统本地时区
    :raises ValueError: 值为空、无法解析、非正数或非有限数字时抛出
    """
    if isinstance(value, bool):
        raise TypeError("bool is not a valid timestamp")
    if isinstance(value, datetime):
        if value.tzinfo is None and tz is not None:
            value = value.replace(tzinfo=tz)
        return int(round(value.timestamp() * NANOSECONDS_PER_SECOND))
    if isinstance(value, (int, float)):
        return _numeric_timestamp_ns(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError("timestamp value must not be an empty string")
        try:
            seconds = float(text)
        except ValueError:
            try:
                dt = datetime.fromisoformat(text.replace("Z", "+00:00"))
            except ValueError as exc:
                raise ValueError(f"Invalid ISO 8601 timestamp: {text!r}") from exc
            if dt.tzinfo is None and tz is not None:
                dt = dt.replace(tzinfo=tz)
            return int(round(dt.timestamp() * NANOSECONDS_PER_SECOND))
        return _numeric_timestamp_ns(seconds)
    raise TypeError(f"Unsupported timestamp type: {type(value).__name__}")


def _numeric_timestamp_ns(value: Union[int, float]) -> int:
    """把秒/毫秒级数字时间戳归一化为纳秒。"""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"Timestamp must be a finite number, got {value!r}")
    # 超过 10 位的视为毫秒（或微秒），逐级归一化到秒
    while abs(value) > TIMESTAMP_S_MAX:
        value = value / 1000
    if value <= 0:
        raise ValueError(f"Expected positive timestamp, got {value!r}")
    return int(round(value * NANOSECONDS_PER_SECOND))
