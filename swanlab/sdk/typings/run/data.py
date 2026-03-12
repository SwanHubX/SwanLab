"""
@author: cunyue
@file: data.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时数据类型
"""

from typing import Literal

MediaTransferType = Literal["image", "audio", "video", "text", "echarts"]


DataTransferType = Literal["scalar", MediaTransferType]
"""
数据类型，包括scalar、image、audio、video、text、table、echarts
"""
