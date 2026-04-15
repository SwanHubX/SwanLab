"""
@author: cunyue
@file: __init__.py
@time: 2026/4/15 23:35
@description: SwanLab SDK cmd 模块的类型定义和导出
"""

from typing import Literal, Union

LoginSaveType = Union[bool, Literal["local", "root"]]
