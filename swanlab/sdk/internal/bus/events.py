"""
@author: cunyue
@file: events.py
@time: 2026/3/13
@description: SwanLab 事件总线协议定义
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, Type, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import UpdateType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.internal.context.transformer import TransformData
from swanlab.sdk.typings.run.column import ScalarXAxisType


@dataclass
class MetricLogEvent:
    """日志记录事件"""

    data: Dict[str, Any]
    step: int
    timestamp: Timestamp


@dataclass
class ScalarDefineEvent:
    """显式创建标量列事件"""

    # 指标键
    key: str
    # 指标名
    name: Optional[str] = None
    # 指标颜色
    color: Optional[str] = None
    # 是否为系统指标
    system: bool = False
    # x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
    x_axis: Optional[ScalarXAxisType] = None
    # 图表索引
    chart: Optional[str] = None
    # 图表名
    chart_name: Optional[str] = None


@dataclass
class ConfigEvent:
    """配置记录事件（path 固定为 files/config.yaml）"""

    update: UpdateType
    timestamp: Timestamp


@dataclass
class ConsoleEvent:
    """控制台输出事件"""

    line: str
    stream: StreamType
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
    ScalarDefineEvent,
    ConfigEvent,
    ConsoleEvent,
    FileSaveEvent,
]

# 数据解析返回类型
ParseResult = Tuple[Optional[Union[MediaRecord, ScalarRecord]], Type[TransformData]]
