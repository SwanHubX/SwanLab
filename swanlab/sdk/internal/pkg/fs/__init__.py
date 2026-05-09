"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 13:36
@description: SwanLab SDK 文件系统辅助函数
这里主要是为了解决NAS等网络文件系统的延迟问题，在python层面这一切被操作系统屏蔽了，因此我们需要一个抗延迟的重试机制
"""

import os
import shutil
from pathlib import Path
from typing import Union

from .dir import safe_mkdir, safe_mkdirs
from .write import safe_write

__all__ = ["safe_mkdir", "safe_mkdirs", "safe_link", "safe_write"]


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
