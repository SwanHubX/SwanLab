"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时类型
"""

from typing import Literal

ResumeType = Literal["must", "allow", "never"]
"""
Run 恢复策略
"""

ModeType = Literal["disabled", "online", "local", "offline"]
"""
Run 运行策略
"""

FinishType = Literal["success", "crashed", "aborted"]
"""
Run 结束策略
"""


AsyncLogType = Literal["asyncio", "threading", "spawn", "fork"]
"""
异步日志上报类型
"""

SaveType = Literal["live", "now", "end"]
"""
文件保存策略
"""
