"""
@author: cunyue
@file: __init__.py
@time: 2026/8/18
@description: 时间戳解析辅助函数，供 API 查询 lane 与运行时 SDK lane 共用
"""

import math
from datetime import datetime, timezone
from decimal import ROUND_HALF_EVEN, Decimal, InvalidOperation
from typing import Union

from google.protobuf.timestamp_pb2 import Timestamp

__all__ = [
    "TimeLike",
    "TIMESTAMP_S_MAX",
    "TIMESTAMP_S_MIN",
    "NANOSECONDS_PER_SECOND",
    "parse_resp_timestamp_s",
    "parse_cli_timestamp_ns",
    "format_time",
    "set_timestamp",
]

# 秒级 Unix 时间戳的合法区间（10 位）
TIMESTAMP_S_MIN = 1_000_000_000
TIMESTAMP_S_MAX = 9_999_999_999

NANOSECONDS_PER_SECOND = 1_000_000_000

_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)

# Unix 时间戳（秒/毫秒）、ISO 8601 字符串或 datetime
TimeLike = Union[int, float, str, datetime]


def parse_resp_timestamp_s(value: Union[int, str]) -> int:
    """解析响应体中的时间戳为秒级 Unix 时间戳。

    支持秒/毫秒级数字时间戳和 ISO 8601 字符串；无时区字符串按 UTC 解释。

    :raises ValueError: 值为空、无法解析、非正数或越界时抛出
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


def parse_cli_timestamp_ns(value: TimeLike) -> int:
    """解析 CLI 输入的时间为纳秒级 Unix 时间戳。

    支持 int / float / 数字字符串 / ISO 8601 字符串 / datetime；
    无时区字符串和 naive datetime 按系统本地时区解释。

    :raises TypeError: 输入类型不支持时抛出
    :raises ValueError: 值为空、无法解析、非有限数字或越界时抛出
    """
    if isinstance(value, bool):
        raise TypeError("bool is not a valid timestamp")
    if isinstance(value, datetime):
        return _datetime_timestamp_ns(value)
    if isinstance(value, (int, float)):
        return _numeric_timestamp_ns(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError("timestamp value must not be an empty string")
        try:
            seconds: Union[int, float] = int(text)
        except ValueError:
            try:
                seconds = float(text)
            except ValueError:
                try:
                    dt = datetime.fromisoformat(text.replace("Z", "+00:00"))
                except ValueError as exc:
                    raise ValueError(f"Invalid ISO 8601 timestamp: {text!r}") from exc
                return _datetime_timestamp_ns(dt)
        return _numeric_timestamp_ns(seconds)
    raise TypeError(f"Unsupported timestamp type: {type(value).__name__}")


def _datetime_timestamp_ns(value: datetime) -> int:
    """把 datetime 转为纳秒时间戳，用整数时间差避免亚秒精度损失。"""
    if value.tzinfo is None:
        # 无时区按本地时区解释
        value = value.astimezone()
    delta = value - _EPOCH
    ns = (delta.days * 86400 + delta.seconds) * NANOSECONDS_PER_SECOND + delta.microseconds * 1000
    if not TIMESTAMP_S_MIN * NANOSECONDS_PER_SECOND <= ns <= TIMESTAMP_S_MAX * NANOSECONDS_PER_SECOND:
        raise ValueError(f"Timestamp {value!r} out of range [{TIMESTAMP_S_MIN}, {TIMESTAMP_S_MAX}] seconds")
    return ns


def _numeric_timestamp_ns(value: Union[int, float]) -> int:
    """把秒/毫秒级数字时间戳归一化为纳秒，用 Decimal 避免浮点误差。"""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"Timestamp must be a finite number, got {value!r}")
    # 超过 10 位的视为毫秒（或微秒），逐级归一化到秒
    try:
        # 用 str 还原浮点字面量
        number = Decimal(str(value)) if isinstance(value, float) else Decimal(value)
    except InvalidOperation as exc:  # pragma: no cover
        raise ValueError(f"Invalid timestamp value: {value!r}") from exc
    while abs(number) > TIMESTAMP_S_MAX:
        number = number.scaleb(-3)
    # 拦截漏写分隔符的日期（如 20261020）
    if not TIMESTAMP_S_MIN <= number <= TIMESTAMP_S_MAX:
        raise ValueError(f"Timestamp {value!r} out of range [{TIMESTAMP_S_MIN}, {TIMESTAMP_S_MAX}] seconds")
    # 已归一化为秒，放大到纳秒
    return int(number.scaleb(9).to_integral_value(rounding=ROUND_HALF_EVEN))


def set_timestamp(target: Timestamp, value: TimeLike) -> None:
    try:
        target.FromNanoseconds(parse_cli_timestamp_ns(value))
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid time value {value!r}: {e}") from e


def format_time(ts: Timestamp) -> str:
    dt = ts.ToDatetime().replace(tzinfo=timezone.utc)
    return dt.astimezone().strftime("%Y-%m-%d %H:%M:%S")
