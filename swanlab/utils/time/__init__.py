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

__all__ = [
    "TIMESTAMP_S_MAX",
    "TIMESTAMP_S_MIN",
    "NANOSECONDS_PER_SECOND",
    "parse_timestamp_s",
    "parse_timestamp_ns",
]

# 秒级 Unix 时间戳的合法区间（10 位）
TIMESTAMP_S_MIN = 1_000_000_000
TIMESTAMP_S_MAX = 9_999_999_999

NANOSECONDS_PER_SECOND = 1_000_000_000

_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)


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


def parse_timestamp_ns(value: Union[int, float, str, datetime]) -> int:
    """将时间统一转换为纳秒级 Unix 时间戳，用于 ``swanlab sync`` 的时间范围过滤。

    支持格式：
    - int / float：秒级或毫秒级 Unix 时间戳，大于 10 位时自动按毫秒归一化为秒
    - 数字字符串：同上
    - ISO 8601 字符串：如 ``"2024-08-01T00:00:00Z"``、``"2024-08-01 08:00:00"``，也支持纯日期
    - ``datetime``：直接取时间戳

    不带时区的字符串和 ``datetime`` 按**系统本地时区**解析，与用户在终端输入墙上时间的习惯一致；
    带时区的输入始终以输入自带的时区为准。

    :param value: 待解析的时间
    :raises ValueError: 值为空、无法解析、非正数或非有限数字时抛出
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
    """把 datetime 转换为纳秒时间戳。

    使用整数时间差而非 ``timestamp() * 1e9``，避免浮点在亚秒位上产生 ns 级误差，
    否则闭区间端点上的记录会被过滤器判到范围外。
    """
    if value.tzinfo is None:
        # 无时区输入按系统本地时区解释，与终端输入习惯一致
        value = value.astimezone()
    delta = value - _EPOCH
    return (delta.days * 86400 + delta.seconds) * NANOSECONDS_PER_SECOND + delta.microseconds * 1000


def _numeric_timestamp_ns(value: Union[int, float]) -> int:
    """把秒/毫秒级数字时间戳归一化为纳秒。

    使用 ``Decimal`` 做整数标度换算，避免 ``value / 1000`` 与 ``* 1e9`` 的浮点误差
    （整数毫秒时间戳经浮点归一化后可能偏出数十到数百纳秒）。
    """
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"Timestamp must be a finite number, got {value!r}")
    # 超过 10 位的视为毫秒（或微秒），逐级归一化到秒
    try:
        # str(float) 还原用户输入的字面量，Decimal(str(x)) 比 Decimal(x) 更接近预期值
        number = Decimal(str(value)) if isinstance(value, float) else Decimal(value)
    except InvalidOperation as exc:  # pragma: no cover - 仅防御非法输入
        raise ValueError(f"Invalid timestamp value: {value!r}") from exc
    while abs(number) > TIMESTAMP_S_MAX:
        number = number.scaleb(-3)
    if number <= 0:
        raise ValueError(f"Expected positive timestamp, got {value!r}")
    # 归一化后 number 已是秒，统一放大到纳秒；Decimal 运算不引入浮点误差
    return int(number.scaleb(9).to_integral_value(rounding=ROUND_HALF_EVEN))
