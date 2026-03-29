"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: SwanLab 事件总线，emitter 与 event 定义独立于 run/，
供 run/ 和 system/ 等模块共同引用，避免循环导入。
"""

from .emitter import RunEmitter, RunQueue
from .events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    EventPayload,
    FlushPayload,
    MetadataEvent,
    MetricLogEvent,
    MetricsUploadEvent,
    MetricsUploadPayload,
    ParseResult,
    RequirementsEvent,
    RunFinishEvent,
    RunStartEvent,
    ScalarDefineEvent,
)

__all__ = [
    "RunEmitter",
    "RunQueue",
    "MetricLogEvent",
    "MetricsUploadEvent",
    "MetricsUploadPayload",
    "ScalarDefineEvent",
    "RunFinishEvent",
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
