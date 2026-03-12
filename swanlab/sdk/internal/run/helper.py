"""
@author: cunyue
@file: helper.py
@time: 2026/3/12 16:38
@description:
"""

import re
from typing import Any, Dict, Mapping, Optional

from swanlab.sdk.internal.pkg import console


def flatten_dict(
    d: Mapping[str, Any], parent_key: str = "", parent_dict: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    展开字典，例如 {"a": {"b": {"c": 1}}} -> {"a/b/c": 1}
    如果出现重复键名，根据顺序，顺序靠后的键名会覆盖靠前的键名
    :param d: 待展开的字典
    :param parent_key: 父级键名，用于构建新键名
    :param parent_dict: 父级字典，用于存储展开后的结果
    :return: 展开后的字典
    """
    # 顶层调用时初始化字典（避免可变默认参数陷阱）
    if parent_dict is None:
        parent_dict = {}

    for k, v in d.items():
        # 防御性编程：用户可能会传非字符串的 key（比如整数），强制转为 str
        k_str = str(k)
        new_key = f"{parent_key}/{k_str}" if parent_key else k_str

        if isinstance(v, Mapping):
            # 递归调用，将同一个 parent_dict 引用传递下去
            flatten_dict(v, new_key, parent_dict)
        else:
            try:
                # 对完整拼接后的路径进行最终的合法性校验、字符替换与截断
                safe_key = validate_key(new_key)
            except Exception as e:
                # 终极容错：如果清洗后变成空字符串（比如用户传了 {"///": 1}），
                # 打印错误并丢弃该字段，绝不中断其他合法数据的解析！
                console.error(f"SwanLab dropped an invalid metric: {e}")
                continue

            # 检查冲突并警告
            if safe_key in parent_dict:
                console.warning(
                    f"Duplicate key found after sanitization: '{safe_key}'. "
                    "The latter value will overwrite the former one."
                )
            # 赋值
            parent_dict[safe_key] = v

    return parent_dict


# 在模块顶层预编译正则，性能更好
_INVALID_KEY_PATTERN = re.compile(r"[^\w./-]")

# 全局警告缓存池，保证相同的非法 key 在同一进程中只警告一次
# 不过在设计上可能每次实验都需要重新初始化更符合直觉一些，不过这属于小概率事件，考虑到性能，将其作为进程级别全局变量
_WARNED_KEYS = set()


def validate_key(key: str, max_len: int = 255) -> str:
    """
    检查并清洗 key 字符串格式。
    将非法字符替换为下划线，自动剥离边缘的非法字符，并在超长时截断。

    :param key: 待检查的键名
    :param max_len: 键名的最大长度，默认为255
    :return: 清洗后的键名
    """
    # 宽容处理类型：如果是 int/float，直接转 str，不抛异常
    if not isinstance(key, str):
        key = str(key)

    original_key = key

    # 剥离头尾的空白字符、'.' 和 '/'
    key = key.strip(" \t\n\r./")

    if not key:
        # 只有在清洗后完全为空这种极端且无法挽救的情况下，才抛出异常
        raise ValueError(f"SwanLab key: '{original_key}' is invalid or empty after sanitization.")

    # 3. 使用预编译的正则进行替换，速度极快
    sanitized_key = _INVALID_KEY_PATTERN.sub("_", key)

    # 长度截断
    if len(sanitized_key) > max_len:
        sanitized_key = sanitized_key[:max_len]

    # 4. 友好提示
    if sanitized_key != original_key:
        if original_key not in _WARNED_KEYS:
            console.warning(
                f"Key '{original_key}' has been sanitized to '{sanitized_key}', please use valid characters."
            )
            _WARNED_KEYS.add(original_key)

    return sanitized_key
