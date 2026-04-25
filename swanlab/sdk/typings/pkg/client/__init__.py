"""
@author: cunyue
@file: __init__.py
@time: 2026/4/14 00:48
@description: SwanLab 客户端模块的类型提示定义
"""

from collections.abc import Mapping, Sequence
from typing import Any, TypeAlias, Union

from .bootstrap import LoginResponse

JSONDict: TypeAlias = Mapping[str, Any]
JSONList: TypeAlias = Sequence[Any]
JSONBody: TypeAlias = Union[JSONDict, JSONList, None]


__all__ = ["JSONBody", "JSONDict", "JSONList", "LoginResponse"]
