"""
@author: cunyue
@file: helper.py
@time: 2026/3/12 16:38
@description:
"""

import re
from typing import Any, Dict, Mapping, Optional, get_args

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.data import ScalarXAxisType


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
    :raises ValueError: 如果清洗后为空字符串
    """
    # 宽容处理类型：如果是 int/float，直接转 str，不抛异常
    if not isinstance(key, str):
        key = str(key)

    original_key = key

    # 剥离头尾的空白字符、'.' 和 '/'
    key = key.strip(" \t\n\r./")

    if not key:
        # 只有在清洗后完全为空这种极端且无法挽救的情况下，才抛出异常
        raise ValueError(
            f"SwanLab key: '{original_key}' is invalid or empty after sanitization, please use valid characters (alphanumeric, '.', '-', '/') and avoid special characters."
        )

    # 3. 使用预编译的正则进行替换，速度极快
    sanitized_key = _INVALID_KEY_PATTERN.sub("_", key)

    # 长度截断
    if len(sanitized_key) > max_len:
        sanitized_key = sanitized_key[:max_len]

    # 4. 友好提示
    if sanitized_key != original_key:
        if original_key not in _WARNED_KEYS:
            console.warning(
                f"Key '{original_key}' has been sanitized to '{sanitized_key}', due to invalid characters or length exceeding limit."
            )
            _WARNED_KEYS.add(original_key)

    return sanitized_key


def safe_validate_log_data(data: Mapping[str, Any]) -> Optional[Mapping[str, Any]]:
    """
    检查并清洗日志数据，如果出现非法键名或值类型不支持，返回 None。

    :param data: 待检查的日志数据
    :return: 清洗后的日志数据或 None
    """
    if not isinstance(data, Mapping):
        return None
    return data


def safe_validate_step(step: Optional[int]) -> Optional[int]:
    """
    检查并清洗 step 整数格式，如果出现非法值，返回 None。

    :param step: 待检查的步数
    :return: 清洗后的步数或 None
    """
    if step is None:
        return None
    if not isinstance(step, int):
        return None
    return step


def safe_validate_key(key: str) -> Optional[str]:
    """
    检查并清洗 key 字符串格式，如果出现非法字符或长度超过限制，返回 None。

    :param key: 待检查的键名
    :return: 清洗后的键名或 None
    """
    try:
        return validate_key(key)
    except ValueError:
        return None


def safe_validate_name(name: Optional[str], max_len: int = 255) -> Optional[str]:
    """
    检查并清洗指标名称，如果出现非法字符或长度超过限制，返回 None。

    :param name: 待检查的指标名称
    :param max_len: 名称的最大长度，默认为255
    :return: 清洗后的指标名称或 None
    """
    if name is None:
        return None
    if not isinstance(name, str):
        return None
    if len(name) > max_len:
        return None
    return name


def safe_validate_chart_name(name: Optional[str], max_len: int = 255) -> Optional[str]:
    """
    检查并清洗图表名称，如果出现非法字符或长度超过限制，返回 None。

    :param name: 待检查的图表名称
    :param max_len: 名称的最大长度，默认为255
    :return: 清洗后的图表名称或 None
    """
    if name is None:
        return None
    if not isinstance(name, str):
        return None
    if len(name) > max_len:
        return None
    return name


def safe_validate_x_axis(x_axis: Optional[ScalarXAxisType]) -> Optional[ScalarXAxisType]:
    """
    检查并清洗 x 轴指标名称，如果出现非法字符或长度超过限制，返回 None。

    :param x_axis: 待检查的 x 轴指标名称
    :return: 清洗后的 x 轴指标名称或 None
    """
    if x_axis is None:
        x_axis = "_step"
    try:
        return validate_key(x_axis)
    except ValueError:
        return None


def safe_validate_color(color: Optional[str]) -> Optional[str]:
    """
    检查并清洗颜色字符串格式，必须是#开头的十六进制颜色代码

    :param color: 待检查的颜色字符串
    :return: 清洗后的颜色字符串或 None
    """
    if color is None:
        return None
    if not color.startswith("#"):
        return None
    if len(color) != 7:
        return None
    if not all(c in "0123456789abcdefABCDEF" for c in color[1:]):
        return None
    return color


def safe_validate_state(state: FinishType) -> Optional[FinishType]:
    """
    检查并清洗运行结束状态，如果出现非法值，返回 None。

    :param state: 待检查的运行结束状态
    :return: 清洗后的运行结束状态或 None
    """
    if state not in get_args(FinishType):
        return None
    return state
