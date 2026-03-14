"""
@author: cunyue
@file: events.py
@time: 2026/3/13
@description: SwanLab 事件总线协议定义
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import UpdateType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.data import DataTransferType, ScalarXAxisType

# ==========================================
# 事件流定义 (Event Bus Definitions)
# ==========================================


@dataclass
class RunStartEvent:
    """运行启动事件"""

    timestamp: Timestamp


@dataclass
class MetricLogEvent:
    """日志记录事件"""

    data: Dict[str, Any]
    step: int
    timestamp: Timestamp


@dataclass
class MetricDefineEvent:
    """显式创建列事件"""

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
    # 指标类型，显式创建目前仅支持标量
    column_type: Literal["scalar"] = "scalar"


@dataclass
class ConfigEvent:
    """配置记录事件（path 固定为 files/config.yaml）"""

    update: UpdateType


@dataclass
class ConsoleEvent:
    """控制台输出事件"""

    line: str
    stream: StreamType


@dataclass
class MetadataEvent:
    """元数据事件"""

    path: str


@dataclass
class RequirementsEvent:
    """依赖记录事件"""

    path: str


@dataclass
class CondaEvent:
    """Conda 环境记录事件"""

    path: str


@dataclass
class RunFinishEvent:
    """运行结束的毒丸信号 (Poison Pill)"""

    state: FinishType
    error: Optional[str]
    timestamp: Optional[Timestamp]


# 事件载体类型
EventPayload = Union[
    MetricLogEvent,
    MetricDefineEvent,
    RunFinishEvent,
    RunStartEvent,
    ConfigEvent,
    ConsoleEvent,
    MetadataEvent,
    RequirementsEvent,
    CondaEvent,
]

# 刷盘载体类型（升级为顶层 Record envelope）
FlushPayload = List[Record]

# 数据解析返回类型
ParseResult = Tuple[Record, DataTransferType]
