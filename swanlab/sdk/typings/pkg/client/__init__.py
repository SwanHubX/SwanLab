"""
@author: cunyue
@file: __init__.py
@time: 2026/4/14 00:48
@description: SwanLab 客户端模块的类型提示定义
"""

import sys
from collections.abc import Mapping, Sequence

if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    from typing_extensions import TypeAlias

from typing import Any, Union

from .bootstrap import LoginResponse

JSONDict: TypeAlias = Mapping[str, Any]
JSONList: TypeAlias = Sequence[Any]
JSONBody: TypeAlias = Union[JSONDict, JSONList, None]


__all__ = ["JSONBody", "JSONDict", "JSONList", "LoginResponse"]
