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


@dataclass
class MetricLogEvent:
    """日志记录事件"""

    data: Dict[str, Any]
    step: int
    timestamp: Timestamp


@dataclass
class MetricDefineEvent:
    """define_metric 事件，携带参数归一化后的定义补丁。

    字段为两态：``None`` 表示 "未提供"，具体值表示使用该值：

    - 公开签名无法区分 "省略" 与 "显式 None"（``x_axis`` / ``section_name``），
      两者均为 ``None``（"未提供"）；merge 时保留旧值，replace 时重置为默认。
      清除旧值必须使用 ``overwrite=True``。
    - ``hidden`` / ``step_sync`` 的 ``False`` 是有效值，不会被收成 "未提供"。
    """

    # rule identity（exact key 或 glob pattern）
    key: str
    # key 是否为末尾单个 '*' 的 glob 模式；由上层校验分类后传入，resolver 不再解析字符串
    is_glob: bool = False
    # None=not provided, str=custom X key
    x_axis: Optional[str] = None
    # None=not provided, str=named section
    section_name: Optional[str] = None
    # None=not provided, True/False
    hidden: Optional[bool] = None
    # None=not provided, True/False。仅对 custom X 有意义：该 Y 是否触发跨 step 的 X 注入
    step_sync: Optional[bool] = None
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
