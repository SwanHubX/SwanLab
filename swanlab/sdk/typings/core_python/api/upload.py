"""
@author: cunyue
@file: upload.py
@time: 2026/4/23 20:01
@description: 上传相关API
"""

from typing import List, Literal, TypedDict


class ConsoleMetric(TypedDict):
    level: Literal["INFO", "ERROR"]
    epoch: int
    message: str
    create_time: str


ConsoleMetrics = List[ConsoleMetric]
"""
Console 日志指标
"""
