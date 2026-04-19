"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: SwanLab 事件总线，emitter 与 event 定义独立于 run/，
供 run/ 和 system/ 等模块共同引用，避免循环导入。
"""

from .emitter import EmitterProtocol, RunEmitter, RunQueue
from .events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    EventPayload,
    FlushPayload,
    MetadataEvent,
    MetricLogEvent,
    ParseResult,
    RequirementsEvent,
    ScalarDefineEvent,
)

__all__ = [
    "RunEmitter",
    "RunQueue",
    "EmitterProtocol",
    "MetricLogEvent",
    "ScalarDefineEvent",
    "ConfigEvent",
    "ConsoleEvent",
    "MetadataEvent",
    "RequirementsEvent",
    "CondaEvent",
    "EventPayload",
    "FlushPayload",
    "ParseResult",
]
