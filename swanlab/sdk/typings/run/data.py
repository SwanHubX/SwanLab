"""
@author: cunyue
@file: data.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时数据类型
"""

from typing import Literal

DataModuleType = Literal["scalar", "image", "audio", "video", "text", "table", "echarts"]
"""
数据类型，包括scalar、image、audio、video、text、table、echarts
"""
