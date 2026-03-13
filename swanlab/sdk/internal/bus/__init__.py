"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: SwanLab 事件总线，emitter 与 event 定义独立于 run/，
供 run/ 和 system/ 等模块共同引用，避免循环导入。
"""

from .emitter import RunEmitter
from .events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    DefineEvent,
    EventPayload,
    FinishEvent,
    FlushPayload,
    LogEvent,
    MetadataEvent,
    ParseResult,
    RequirementsEvent,
    RunStartEvent,
)

__all__ = [
    "RunEmitter",
    "LogEvent",
    "DefineEvent",
    "FinishEvent",
    "RunStartEvent",
    "ConfigEvent",
    "ConsoleEvent",
    "MetadataEvent",
    "RequirementsEvent",
    "CondaEvent",
    "EventPayload",
    "FlushPayload",
    "ParseResult",
]
