"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 19:08
@description: SwanLab SDK 格式化工具
"""

import re
from typing import List


def validate_project(name: str, max_len: int = 100) -> str:
    """
    校验项目名称
    """
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Project name must be a non-empty string.")

    name = name.strip()
    if len(name) > max_len:
        raise ValueError(f"Project name '{name}' exceeds the maximum length of {max_len} characters.")

    if not re.match(r"^[0-9a-zA-Z_\-+.]+$", name):
        raise ValueError(f"Project name '{name}' contains invalid characters. Allowed: 0-9, a-z, A-Z, _, -, +, .")

    return name


def validate_experiment(name: str, max_len: int = 250) -> str:
    """
    校验实验名称
    """
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Experiment name must be a non-empty string.")

    name = name.strip()
    if len(name) > max_len:
        raise ValueError(f"Experiment name '{name}' exceeds the maximum length of {max_len} characters.")

    return name


def validate_description(desc: str, max_len: int = 1024) -> str:
    """
    校验描述信息
    """
    if not isinstance(desc, str):
        raise TypeError("Description must be a string.")

    if len(desc) > max_len:
        raise ValueError(f"Description exceeds the maximum length of {max_len} characters.")

    return desc


def validate_tags(tags: List[str], max_len: int = 200) -> List[str]:
    """
    校验标签列表
    """
    if not isinstance(tags, list):
        raise TypeError("Tags must be a list of strings.")

    for tag in tags:
        if not isinstance(tag, str):
            raise TypeError(f"Tag '{tag}' must be a string.")
        if len(tag) > max_len:
            raise ValueError(f"Tag '{tag}' exceeds the maximum length of {max_len} characters.")

    return tags


def validate_job_type(job_type: str, max_len: int = 256) -> str:
    """
    校验任务类型
    """
    if not isinstance(job_type, str):
        raise TypeError("Job type must be a string.")

    if len(job_type) > max_len:
        raise ValueError(f"Job type '{job_type}' exceeds the maximum length of {max_len} characters.")

    return job_type


def validate_group(group: str, max_len: int = 256) -> str:
    """
    校验实验组
    """
    if not isinstance(group, str):
        raise TypeError("Group must be a string.")

    if len(group) > max_len:
        raise ValueError(f"Group '{group}' exceeds the maximum length of {max_len} characters.")

    return group


def validate_run_id(run_id: str, max_len: int = 64) -> str:
    """
    校验运行 ID
    """
    if not isinstance(run_id, str):
        raise TypeError(f"Run ID must be a string, got {type(run_id).__name__}.")

    run_id = run_id.strip()

    # 使用 f-string 动态构建正则。注意大括号 {{ 和 }} 是用来在 f-string 中转义普通大括号的
    pattern = rf"^[a-z0-9_]{{1,{max_len}}}$"
    if not re.match(pattern, run_id):
        raise ValueError(
            f"Run ID '{run_id}' is invalid. "
            f"It must be between 1 and {max_len} characters long, containing only lowercase letters, digits, and underscores (_)."
        )

    return run_id


def validate_metric_key(key: str, max_len: int = 255) -> str:
    """
    校验指标/日志键名 (Key)：
    默认最大长度 255 个字符，不能为空，且不能以 `.` 或 `/` 开头或结尾。
    :param key: 指标/日志键名
    :param max_len: 最大允许的字符串长度
    """
    if not isinstance(key, str):
        raise TypeError(f"Key '{key}' is not a string.")

    key = key.strip()
    if not key:
        raise ValueError("Key cannot be an empty string.")

    if len(key) > max_len:
        raise ValueError(f"Key '{key}' exceeds the maximum length of {max_len} characters.")

    if key.startswith((".", "/")) or key.endswith((".", "/")):
        raise ValueError(f"Key '{key}' cannot start or end with '.' or '/'.")

    return key
