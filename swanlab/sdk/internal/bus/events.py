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

    字段使用 ``UNSET`` 哨兵区分 "未提供" 与显式值：

    - ``UNSET``: 参数未提供。公开 API 的省略与显式 ``None`` / ``False``
      在归一化阶段统一转为 ``UNSET``（公开签名无法区分 "省略" 与
      "显式 None"）；merge 时保留旧值，replace 时重置为默认。
    - ``x_axis`` / ``section_name``: ``UNSET`` 或具体值。清除旧值必须
      使用 ``overwrite=True``；resolver 对 ``None`` 的处理分支当前不可达，
      仅为数据模型完备性保留。
    - ``hidden``: ``UNSET`` 或 ``True``（公开层 ``False`` 亦收成 ``UNSET``）。
    - ``step_sync``: ``UNSET`` / ``True`` / ``False``。``False`` 是有效值，
      不可收成 ``UNSET``。
    """

    # rule identity（exact key 或 glob pattern）
    key: str
    # UNSET=未提供（含公开层显式 None），str=custom X key；None 分支当前不可达
    x_axis: Any = UNSET
    # UNSET=未提供（含公开层显式 None），str=named section；None 分支当前不可达
    section_name: Any = UNSET
    # UNSET=未提供（公开层 False 亦收成 UNSET），True=隐藏
    hidden: Any = UNSET
    # UNSET=未提供，True/False。仅对 custom X 有意义：该 Y 是否触发跨 step 的 X 注入
    step_sync: Any = UNSET
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
