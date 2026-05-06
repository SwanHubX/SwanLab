"""
@author: cunyue
@file: system.py
@time: 2026/5/5 19:31
@description: 系统指标
"""

__all__ = ["fmt_system_key", "is_system_key"]

_SYSTEM_KEY_PREFIX = "__swanlab__."


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
