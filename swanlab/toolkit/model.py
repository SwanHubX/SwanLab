"""
@author: cunyue
@file: model.py
@time: 2025/7/20 16:25
@description: 部分新增的模型定义（patch)
"""

from typing import TypedDict


class LogContent(TypedDict):
    """日志内容字典类型

    结构示例:
    {
        "message": "hello world",
        "create_time": "2025-05-15 18:35:00",
        "epoch": 1
    }
    """

    message: str
    create_time: str
    epoch: int
