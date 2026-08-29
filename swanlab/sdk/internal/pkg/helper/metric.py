"""
@author: cunyue
@file: metric.py
@time: 2026/5/5 19:31
@description: 指标 key 工具（系统指标 key 前缀、系统 X 轴的判定）
"""

__all__ = ["fmt_system_key", "is_system_key", "is_custom_x"]

_SYSTEM_KEY_PREFIX = "__swanlab__."

# 系统 X 轴：由服务端/SDK 内部维护，不视为自定义 metric key
_SYSTEM_X_AXES = frozenset({"_step", "_relative_time"})


def fmt_system_key(key: str):
    """
    格式化系统指标key，如果已经是系统指标key，报错
    """
    if key.startswith(_SYSTEM_KEY_PREFIX):
        raise ValueError(f"System metric key '{key}' is already a system metric key")
    return f"{_SYSTEM_KEY_PREFIX}{key}"


def is_system_key(key: str):
    """
    判断是否是系统指标key
    :param key: 指标key
    :return: 是否是系统指标key
    """
    return key.startswith(_SYSTEM_KEY_PREFIX)


def is_custom_x(x_axis: str) -> bool:
    """
    判断 x_axis 是否为自定义（非系统 _step / _relative_time）的 scalar key
    :param x_axis: X 轴 key
    :return: 是否为自定义 X 轴
    """
    return x_axis not in _SYSTEM_X_AXES
