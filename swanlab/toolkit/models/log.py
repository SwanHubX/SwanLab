"""
@author: cunyue
@file: log.py
@time: 2025/10/10 22:01
@description: 终端日志模型
"""

from typing import TypedDict

__all__ = ['LogContent']


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
