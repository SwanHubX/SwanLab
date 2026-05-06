"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:18
@description: SwanLab SDK 辅助函数
"""

from pathlib import PurePosixPath

from .env import DEBUG, is_interactive, is_jupyter
from .system import fmt_system_key, is_system_key
from .version import get_swanlab_latest_version, get_swanlab_version

__all__ = [
    "DEBUG",
    "is_jupyter",
    "is_interactive",
    "strip_none",
    "fmt_run_path",
    "fmt_system_key",
    "is_system_key",
    "get_swanlab_version",
    "get_swanlab_latest_version",
]


def strip_none(data: dict, strip_empty_dict: bool = True, strip_empty_str: bool = False) -> dict:
    """
    递归剔除字典中的 None 值。

    :param data: 待处理的字典
    :param strip_empty_dict: 是否连同空字典一起剔除（包含剔除 None 后变为空的嵌套字典），默认为 True
    :param strip_empty_str: 是否连同空字符串一起剔除，仅当 strip_empty_dict 为 True 时有效，默认为 False
    """
    clean_data = {}
    for k, v in data.items():
        if isinstance(v, dict):
            # 递归调用时，记得把参数透传下去
            cleaned_v = strip_none(v, strip_empty_dict=strip_empty_dict, strip_empty_str=strip_empty_str)

            # 只有当嵌套字典里有值，或者我们明确声明“不剔除空字典”时，才保留这个 key
            if cleaned_v or not strip_empty_dict:
                clean_data[k] = cleaned_v
        elif v is not None and (not strip_empty_str or v != ""):
            clean_data[k] = v

    return clean_data


def fmt_run_path(run_path: str) -> str:
    """
    格式化运行路径为 /@:username/:project_name/runs/:run_id，适配前端路由设计
    传入路径格式为 /:username/:project_name/:run_id
    后端返回可能是：/cunyue/demo/abc123，cunyue/demo/abc123，/@cunyue/demo/abc123，@cunyue/demo/abc123
    """
    parts = PurePosixPath(run_path).parts

    # 去掉根路径 "/"
    parts = parts[1:] if parts and parts[0] == "/" else parts

    if len(parts) < 3:
        raise ValueError(f"Invalid run path: {run_path!r}")

    username, project_name, run_id = parts[0], parts[1], parts[2]

    if not username.startswith("@"):
        username = f"@{username}"

    return PurePosixPath("/", username, project_name, "runs", run_id).as_posix()
