"""
@author: cunyue
@file: __init__.py
@time: 2026/7/29
@description: mimetype 封装，注册自定义扩展并提供类型推断
"""

from __future__ import annotations

import functools
import os
from typing import Union

__all__ = ["guess_type"]

_DEFAULT_MIME_TYPE = "application/octet-stream"


@functools.cache
def _ensure_mimetypes():
    """惰性导入 mimetypes 并注册自定义扩展，仅首次调用时执行。"""
    import mimetypes as _std_mimetypes

    # 注册扩展到标准类型表（guess_type 默认只查标准表）
    # strict=True 会覆盖系统已有映射（源码直接赋值，不抛错），保证跨平台确定性
    _std_mimetypes.add_type("application/yaml", ".yaml")
    _std_mimetypes.add_type("application/yaml", ".yml")

    # 图像类型：与前端 image/* 渲染器对齐，避免各平台 mime 差异
    _std_mimetypes.add_type("image/png", ".png")
    _std_mimetypes.add_type("image/jpeg", ".jpg")
    _std_mimetypes.add_type("image/jpeg", ".jpeg")
    _std_mimetypes.add_type("image/gif", ".gif")
    _std_mimetypes.add_type("image/bmp", ".bmp")
    _std_mimetypes.add_type("image/webp", ".webp")
    _std_mimetypes.add_type("image/svg+xml", ".svg")
    _std_mimetypes.add_type("image/svg+xml", ".svgz")
    _std_mimetypes.add_type("image/tiff", ".tiff")
    _std_mimetypes.add_type("image/tiff", ".tif")

    # 代码 / 文本类型：与前端 CodeRender / MarkdownRender 对齐
    _std_mimetypes.add_type("text/x-python", ".py")
    _std_mimetypes.add_type("text/plain", ".txt")
    _std_mimetypes.add_type("application/json", ".json")
    _std_mimetypes.add_type("text/markdown", ".md")
    _std_mimetypes.add_type("text/markdown", ".markdown")

    return _std_mimetypes


def guess_type(path: Union[str, "os.PathLike[str]"]) -> str:
    """根据文件路径推断 MIME 类型，无法识别时回退到 application/octet-stream。"""
    return _ensure_mimetypes().guess_type(os.fspath(path))[0] or _DEFAULT_MIME_TYPE
