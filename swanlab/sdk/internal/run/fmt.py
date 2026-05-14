"""
@author: cunyue
@file: utils_fmt.py
@time: 2026/3/12 16:38
@description: 一些格式化工具
"""

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple, Union, get_args

from pydantic import ValidationError

from swanlab.sdk.internal.pkg import console, constraints, safe
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.column import ScalarXAxisType


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
            # 如果清洗后变成空字符串（比如用户传了 {"///": 1}），丢弃这个字段
            with safe.block(message="SwanLab dropped an invalid metric"):
                # 对完整拼接后的路径进行最终的合法性校验、字符替换与截断
                safe_key = validate_key(new_key)
                # 检查冲突并警告
                if safe_key in parent_dict:
                    console.warning(
                        f"Duplicate key found after sanitization: '{safe_key}'. "
                        "The latter value will overwrite the former one."
                    )
                # 赋值
                parent_dict[safe_key] = v

    return parent_dict


# 全局警告缓存池，保证相同的非法 key 在同一进程中只警告一次
# 不过在设计上可能每次实验都需要重新初始化更符合直觉一些，不过这属于小概率事件，考虑到性能，将其作为进程级别全局变量
_WARNED_KEYS = set()


def validate_key(key: str, max_len: int = constraints.METRIC_KEY_MAX_LENGTH) -> str:
    """
    检查并清洗 key 字符串格式。
    将非法字符替换为下划线，自动剥离边缘的非法字符，并在超长时截断。

    :param key: 待检查的键名
    :param max_len: 键名的最大长度，默认为255
    :return: 清洗后的键名
    :raises ValueError: 如果清洗后为空字符串或者无法通过 MetricKey 校验
    """
    # 宽容处理类型：如果是 int/float，直接转 str，不抛异常
    if not isinstance(key, str):
        key = str(key)

    max_len = min(max_len, constraints.METRIC_KEY_MAX_LENGTH)

    original_key = key

    # 剥离头尾的空白字符、'.' 和 '/'
    key = key.strip(" \t\n\r./")

    if not key:
        # 只有在清洗后完全为空这种极端且无法挽救的情况下，才抛出异常
        raise ValueError(
            f"SwanLab key: '{original_key}' is invalid or empty after sanitization, please use valid characters and avoid leading/trailing special characters."
        )

    # 长度截断
    if len(key) > max_len:
        key = key[:max_len]

    # 友好提示（仅因头尾剥离或截断而改变）
    if key != original_key:
        if original_key not in _WARNED_KEYS:
            console.warning(
                f"Key '{original_key}' has been trimmed to '{key}', due to leading/trailing characters or length exceeding limit."
            )
            _WARNED_KEYS.add(original_key)

    # 使用 MetricKey 约束做最终校验；内部含控制字符等非法内容时直接抛出，不做替换
    try:
        return constraints.ta_metric_key.validate_python(key)
    except ValidationError as e:
        raise ValueError(str(e)) from e


def safe_validate_log_data(data: Mapping[str, Any]) -> Optional[Mapping[str, Any]]:
    """
    检查并清洗日志数据，如果出现非法键名或值类型不支持，返回 None。

    :param data: 待检查的日志数据
    :return: 清洗后的日志数据或 None
    """
    if not isinstance(data, Mapping):
        return None
    return data


def safe_validate_key(key: str) -> Optional[str]:
    """
    检查 key 字符串格式，如果出现非法字符或长度超过限制，返回 None。

    :param key: 待检查的键名
    :return: 合法的键名或 None
    """
    try:
        return constraints.ta_metric_key.validate_python(key)
    except ValidationError:
        return None


def safe_validate_name(name: Optional[str]) -> Optional[str]:
    """
    检查并清洗指标名称，如果出现非法字符或长度超过限制，返回 None。

    :param name: 待检查的指标名称
    :return: 清洗后的指标名称或 None
    """
    if name is None:
        return None
    try:
        return constraints.ta_metric_name.validate_python(name)
    except ValidationError:
        return None


def safe_validate_chart_name(name: Optional[str]) -> Optional[str]:
    """
    检查并清洗图表名称，如果出现非法字符或长度超过限制，返回 None。

    :param name: 待检查的图表名称
    :return: 清洗后的图表名称或 None
    """
    if name is None:
        return None
    try:
        return constraints.ta_chart_name.validate_python(name)
    except ValidationError:
        return None


def safe_validate_x_axis(x_axis: Optional[ScalarXAxisType]) -> Optional[ScalarXAxisType]:
    """
    检查并清洗 x 轴指标名称，如果出现非法字符或长度超过限制，返回 None。

    :param x_axis: 待检查的 x 轴指标名称
    :return: 清洗后的 x 轴指标名称或 None
    """
    if x_axis is None:
        x_axis = "_step"
    return safe_validate_key(x_axis)


def safe_validate_color(color: Optional[str]) -> Optional[str]:
    """
    检查并清洗颜色字符串格式，必须是#开头的十六进制颜色代码

    :param color: 待检查的颜色字符串
    :return: 清洗后的颜色字符串或 None
    """
    if color is None:
        return None
    try:
        return constraints.ta_hex_color.validate_python(color)
    except ValidationError:
        return None


def safe_validate_state(state: FinishType) -> Optional[FinishType]:
    """
    检查并清洗运行结束状态，如果出现非法值，返回 None。

    :param state: 待检查的运行结束状态
    :return: 清洗后的运行结束状态或 None
    """
    if state not in get_args(FinishType):
        return None
    return state


def _infer_save_base_path(glob_path: Path) -> Path:
    """从绝对 glob 路径推断默认 base_path。

    沿路径各段向前扫描，遇到首个包含 glob 通配符 (``*?["]``) 的段时停止，
    取稳定段的父目录作为 base_path，以保留最近的稳定目录层级。

    例如 ``/mnt/folder/**/*.pt`` → base 为 ``/mnt``，
    相对路径保留 ``folder/**/*.pt``。

    :param glob_path: 已 resolve 的绝对 glob 路径。
    :return: 推断出的 base_path。
    """
    stable_parts: list[str] = []
    for part in glob_path.parts:
        if any(c in part for c in "*?["):
            break
        stable_parts.append(part)

    if len(stable_parts) == len(glob_path.parts):
        anchor = glob_path.parent
    else:
        anchor = Path(*stable_parts)
    return anchor.parent


def resolve_save_paths(
    glob_str: Union[str, bytes],
    base_path: Optional[Union[str, Path]] = None,
) -> Optional[Tuple[Path, Path]]:
    """Validate and resolve glob and base paths for swanlab.save().

    Resolves the glob pattern and base path against the real filesystem,
    validates that the resolved glob does not escape the base, and returns
    the resolved pair. Returns ``None`` if the input is invalid.

    :param glob_str: Glob pattern matching files to save.
    :param base_path: Base directory for relative path resolution.
        Defaults to cwd when the pattern is relative; for absolute patterns
        the parent directory is used and a warning is printed.
    :return: ``(resolved_glob, resolved_base)`` or ``None``.
    """
    if isinstance(glob_str, bytes):
        glob_str = glob_str.decode("utf-8")

    if glob_str.startswith(("gs://", "s3://")):
        console.warning(f"'{glob_str}' is a cloud storage URL and cannot be saved to SwanLab.")
        return None

    glob_path = Path(glob_str)
    resolved_glob = glob_path.resolve()

    if base_path is not None:
        resolved_base = Path(base_path).resolve()
    elif not glob_path.is_absolute():
        resolved_base = Path(".").resolve()
    else:
        resolved_base = _infer_save_base_path(resolved_glob)
        console.warning(
            f"Saving files from absolute path '{glob_str}' without a base_path. "
            f"SwanLab will default to base_path='{resolved_base}' "
            "to preserve the immediate parent directory. "
            "If you want to preserve a different directory structure, "
            "explicitly pass 'base_path' to swanlab.save"
        )

    try:
        resolved_glob.relative_to(resolved_base)
    except ValueError:
        console.error(
            f"Glob pattern '{glob_str}' resolves to '{resolved_glob}', "
            f"which is outside the base path '{resolved_base}'."
        )
        return None

    return resolved_glob, resolved_base
