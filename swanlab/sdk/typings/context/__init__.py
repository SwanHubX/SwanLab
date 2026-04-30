"""
@author: cunyue
@file: __init__.py
@time: 2026/4/30
@description: SwanLab SDK context 模块的类型定义和导出
"""

from typing import Iterable, Union

from swanlab.sdk.protocol import Callback

CallbacksType = Union[Iterable[Callback], Callback]
"""
swanlab 回调函数参数类型，支持单个 Callback 对象或可迭代对象
"""
