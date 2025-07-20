"""
@author: cunyue
@file: types.py
@time: 2025/5/17 17:31
@description: 定义相关类型
"""

from typing import TypedDict, Literal, Callable, List, Any

from swanlab.toolkit import LogContent

# 支持的代理类型
ProxyType = Literal['all', 'stdout', 'stderr', 'none']
# 支持的日志类型
LogType = Literal['stdout', 'stderr']


class LogData(TypedDict):
    """日志数据字典类型

    结构示例:
    {
        "type": "stdout",  # 或 "stderr"
        "content": [{
            "message": "hello world",
            "create_time": "2025-05-15 18:35:00",
            "epoch": 1
        }]
    }
    """

    type: LogType
    contents: List[LogContent]


# 日志写入器类型
WriteHandler = Callable[[str], None]
# 日志处理器类型
LogHandler = Callable[[LogData], Any]
