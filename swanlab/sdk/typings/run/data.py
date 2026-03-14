"""
@author: cunyue
@file: data.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时数据类型
"""

from typing import Literal, Union

MediaTransferType = Literal["image", "audio", "video", "text", "echarts"]
"""
媒体数据类型
"""


DataTransferType = Literal["scalar", MediaTransferType]
"""
数据类型，包括媒体数据类型
"""

ScalarXAxisType = Union[Literal["_step", "_relative_time"], str]
