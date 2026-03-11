"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 13:06
@description: SwanLab 数据转换模块，将用户输入的数据封装为Protobuf格式
"""

from .scalar import Scalar
from .text import Text

__all__ = ["Text", "Scalar"]
