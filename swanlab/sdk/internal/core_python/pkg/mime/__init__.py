"""
@author: cunyue
@file: __init__.py
@time: 2026/7/29
@description: mimetype 封装，注册自定义扩展并提供类型推断
"""

from __future__ import annotations

import mimetypes as _std_mimetypes
import os
from typing import Union

__all__ = ["guess_type"]

_DEFAULT_MIME_TYPE = "application/octet-stream"

# 注册 YAML 扩展到标准类型表（guess_type 默认只查标准表）
# strict=True 会覆盖系统已有映射（源码直接赋值，不抛错），保证跨平台确定性
_std_mimetypes.add_type("application/yaml", ".yaml")
_std_mimetypes.add_type("application/yaml", ".yml")


def guess_type(path: Union[str, "os.PathLike[str]"]) -> str:
    """根据文件路径推断 MIME 类型，无法识别时回退到 application/octet-stream。"""
    return _std_mimetypes.guess_type(os.fspath(path))[0] or _DEFAULT_MIME_TYPE
