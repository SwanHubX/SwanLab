"""
@author: cunyue
@file: events.py
@time: 2026/3/13
@description: SwanLab 事件总线协议定义
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Type, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel
from swanlab.sdk.internal.context.transformer import TransformData

# ── define_metric 用的 presence 哨兵 ──────────────────────────────────────────


class _UnsetType:
    """Sentinel for "field not provided" in definition patches.

    Distinct from ``None`` (which means "clear to system default")
    and ``False``. Used by :class:`MetricDefineEvent` and the resolver's
    :class:`DefinitionPatch` to implement merge vs. replace semantics.
    """

    _instance: Optional["_UnsetType"] = None

    def __new__(cls) -> "_UnsetType":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:
        return "UNSET"

    def __bool__(self) -> bool:
        return False


UNSET: Any = _UnsetType()


@dataclass
class MetricLogEvent:
    """日志记录事件"""

    data: Dict[str, Any]
    step: int
    timestamp: Timestamp


@dataclass
class MetricDefineEvent:
    """define_metric 事件，携带参数归一化后的定义补丁。

    字段使用 ``UNSET`` 哨兵区分 "未提供" 与显式 ``None`` / ``False``：

    - ``UNSET``: 用户未提供该参数；merge 时保留旧值，replace 时重置为默认。
    - ``None`` (x_axis / section_name): 显式清除为系统默认。
    - 具体值: 使用该值。

    注意：公开 API 的 ``x_axis=None`` / ``section_name=None`` 在归一化阶段统一转为 ``UNSET``，
    因为公开签名无法区分 "省略" 与 "显式 None"。清除旧值必须使用 ``overwrite=True``。
    """

    # rule identity（exact key 或 glob pattern）
    key: str
    # UNSET=not provided, None=system default (_step), str=custom X key
    x_axis: Any = UNSET
    # UNSET=not provided, None=default section, str=named section
    section_name: Any = UNSET
    # UNSET=not provided, True/False
    hidden: Any = UNSET
    # merge（False）vs replace（True）
    overwrite: bool = False


@dataclass
class ConfigEvent:
    """配置记录事件"""

    path: Path
    timestamp: Timestamp


@dataclass
class LogEvent:
    """控制台输出事件"""

    line: str
    level: LogLevel
    timestamp: Timestamp


@dataclass
class FileSaveEvent:
    """文件保存事件"""

    source_path: str
    name: str
    policy: str  # "now" | "end" | "live"


# 事件载体类型
EventPayload = Union[
    MetricLogEvent,
    MetricDefineEvent,
    ConfigEvent,
    LogEvent,
    FileSaveEvent,
]

# 数据解析返回类型
ParseResult = Tuple[Optional[Union[MediaRecord, ScalarRecord]], Type[TransformData]]
