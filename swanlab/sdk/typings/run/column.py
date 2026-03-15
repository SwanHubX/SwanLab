"""
@author: cunyue
@file: data.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时数据类型
"""

from typing import Literal, Union

ScalarXAxisType = Union[Literal["_step", "_relative_time"], str]
