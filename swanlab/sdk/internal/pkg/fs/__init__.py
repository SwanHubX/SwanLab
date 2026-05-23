"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 13:36
@description: SwanLab SDK 文件系统辅助函数
这里主要是为了解决NAS等网络文件系统的延迟问题，在python层面这一切被操作系统屏蔽了，因此我们需要一个抗延迟的重试机制
"""

import os
import re
import shutil
from pathlib import Path
from typing import Tuple, Union

from .. import constraints
from .dir import safe_mkdir, safe_mkdirs
from .write import safe_write

__all__ = ["safe_fmt", "safe_mkdir", "safe_mkdirs", "safe_link", "safe_truncate", "safe_write"]


def safe_link(source: Union[str, Path], target: Union[str, Path]) -> Path:
    """将 source 链接到 target（优先 symlink，失败则 copy2 fallback）。自动创建目标父目录。"""
    source, target = Path(source), Path(target)
    safe_mkdir(target.parent)
    if not target.exists():
        try:
            os.symlink(str(source), str(target))
        except OSError:
            shutil.copy2(str(source), str(target))
    return target


_OS_SAFE_RE = re.compile(constraints.OS_SAFE_PATTERN)


def safe_fmt(dirname: str, fallback: str = "unknown") -> str:
    """
    格式化目录名，将不安全字符、"." 以及 "-" 替换为 "_"

    :param dirname: 原始目录名
    :param fallback: 格式化后为空时的兜底值
    :return: 格式化后的目录名
    """
    dirname = "".join(c if _OS_SAFE_RE.match(c) else "_" for c in dirname)
    dirname = dirname.replace(".", "_").replace("-", "_")
    dirname = re.sub(r"_{2,}", "_", dirname)
    dirname = dirname.strip().strip("_")
    return dirname or fallback


def safe_truncate(name: str, max_length: int = 255) -> Tuple[str, bool]:
    """
    截断目录名，保留前3和后部分字符，中间用 "..." 连接，确保结果的 UTF-8 字节长度不超过 max_length。
    例如： "abcdefghijklmnopqrstuvwxyz" 截断为 "abc...yz"（max_length >= 8 时）

    :param name: 原始目录名
    :param max_length: 目录名最大字节长度，必须 >= 7 以容纳 "abc...yz" 格式，一般文件系统的限制是不超过255字节
    :return: 截断后的目录名和是否发生了截断
    """
    if len(name.encode("utf-8")) <= max_length:
        return name, False
    if max_length < 7:
        raise ValueError(f"max_length ({max_length}) is too small for truncate_dirname, must be >= 7.")
    tail_len = max_length - 6
    truncated_name = f"{name[:3]}...{name[-tail_len:]}"
    while len(truncated_name.encode("utf-8")) > max_length and tail_len > 0:
        tail_len -= 1
        truncated_name = f"{name[:3]}...{name[-tail_len:]}" if tail_len > 0 else f"{name[:3]}..."
    return truncated_name, True
