"""
@author: cunyue
@file: __init__.py
@time: 2026/3/12
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理 (基于 Event-Bus 事件驱动架构)
2. 运行、实验上下文管理
3. 触发异步微批处理落盘与回调
"""

import queue
import threading
from dataclasses import dataclass
from functools import cached_property, singledispatchmethod
from pathlib import Path
from typing import Any, Dict, List, Literal, Mapping, Optional, Tuple, Union, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.data.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.data.v1.log_pb2 import LogRecord
from swanlab.proto.swanlab.data.v1.scalar_pb2 import ScalarValue
from swanlab.sdk.internal.context import RunContext, RunMetrics
from swanlab.sdk.internal.context.callbacker import CallbackManager
from swanlab.sdk.internal.context.transformer import TransformMediaType
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType

from .callbackers import CloudCallback, LocalCallback, OfflineCallback
from .data.transforms import Scalar, Text
from .helper import (
    flatten_dict,
    safe_validate_chart_name,
    safe_validate_color,
    safe_validate_key,
    safe_validate_name,
    safe_validate_x_axis,
)

__all__ = ["SwanLabRun", "CloudCallback", "LocalCallback", "OfflineCallback"]

from ...typings.run.data import DataTransferType, ScalarXAxisType
from ..pkg.fs import safe_mkdir

# ==========================================
# 事件流定义 (Event Bus Definitions)
# ==========================================


@dataclass
class LogEvent:
    """日志记录事件"""

    data: Dict[str, Any]
    step: int
    timestamp: Timestamp


@dataclass
class DefineEvent:
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
class FinishEvent:
    """运行结束的毒丸信号 (Poison Pill)"""

    pass


# 事件载体类型
EventPayload = Union[LogEvent, DefineEvent, FinishEvent]

# 刷盘载体类型
FlushPayload = List[Union[LogRecord, ColumnRecord]]

# 数据解析返回类型
ParseResult = Tuple[LogRecord, DataTransferType]


class SwanLabRun:
    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx

        # 事件总线：最大积压 10 万个事件（防 OOM 兜底）
        self._queue: queue.Queue[EventPayload] = queue.Queue(maxsize=100_000)

        # 启动唯一的后台消费者线程（非守护线程，保证安全落盘）
        self._background_logger = threading.Thread(
            target=self._background_log, name="SwanLab-Background-Logger", daemon=True
        )
        self._background_logger.start()

    # ----------------------------------
    # 属性 (Properties)
    # ----------------------------------

    @cached_property
    def id(self) -> str:
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def run_dir(self) -> Path:
        assert self._ctx.run_dir is not None, "Run dir is not set."
        return self._ctx.run_dir

    @cached_property
    def _flush_file(self) -> Path:
        assert self.run_dir is not None, "Run dir is not set when attempting to access _flush_file."
        return self._ctx.backup_file

    @cached_property
    def _callbacker(self) -> CallbackManager:
        return self._ctx.callbacker

    @cached_property
    def _metrics(self) -> RunMetrics:
        return self._ctx.metrics

    # ----------------------------------
    # 生产者 API：只负责发事件，绝不阻塞主线程业务逻辑
    # ----------------------------------

    def log(self, data: Mapping[str, Any], step: Optional[int] = None):
        """记录一组日志（可能触发隐式列创建）"""
        if not isinstance(data, Mapping):
            console.error(f"Log data must be a dict, but got {type(data).__name__}. SwanLab will ignore it.")
            return
        if not isinstance(step, int) and step is not None:
            console.error(f"Step must be an integer or None, but got {type(step).__name__}. SwanLab will ignore it.")
            return

        step = self._metrics.next_step(step)

        ts = Timestamp()
        ts.GetCurrentTime()

        # 展平字典并在内部进行合规性验证和截断
        safe_data = flatten_dict(data)

        # 推送日志事件
        self._queue.put(LogEvent(data=safe_data, step=step, timestamp=ts), block=True)

    def log_text(self, key: str, data: Union[str, Text], caption: Optional[str] = None, step: Optional[int] = None):
        """
        A syntactic sugar for logging text data.
        :param key: The key for the text data.
        :param data: The text data itself or a Text object.
        :param caption: Optional caption for the text data.
        :param step: Optional step for the text data.
        """
        if not isinstance(data, Text):
            data = Text(data, caption=caption)
        self.log({key: data}, step=step)

    def define_scalar(
        self,
        key: str,
        name: Optional[str] = None,
        color: Optional[str] = None,
        x_axis: Optional[ScalarXAxisType] = None,
        chart_name: Optional[str] = None,
    ):
        """
        Explicitly define a scalar column.
        :param key: The key for the scalar column.
        :param name: Optional name for the scalar column.
        :param color: Optional color for the scalar column.
        :param x_axis: Optional x-axis for the scalar column.
        :param chart_name: Optional name for the chart.
        """
        if not (this_key := safe_validate_key(key)):
            return console.error(
                f"Invalid key for define scalar: {key}, please use valid characters (alphanumeric, '.', '-', '/') and avoid special characters."
            )

        original_name = name
        if name and not (name := safe_validate_name(name)):
            return console.error(f"Invalid name for define scalar: {original_name}, must be a string.")

        original_color = color
        if color and not (color := safe_validate_color(color)):
            return console.error(f"Invalid color for define scalar: {original_color}, must be a hex color code.")

        if (this_x_axis := safe_validate_x_axis(x_axis)) is None:
            return console.error(f"Invalid x_axis for define scalar: {x_axis}, must be a valid ScalarXAxisType.")

        original_chart_name = chart_name
        if chart_name and not (chart_name := safe_validate_chart_name(chart_name)):
            return console.error(f"Invalid chart_name for define scalar: {original_chart_name}, must be a string.")

        self._define_scalar(
            key=this_key,
            name=name,
            color=color,
            x_axis=this_x_axis,
            system=False,
            chart_name=chart_name,
        )

    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """安全关闭当前 Run，等待所有日志落盘"""
        state = state.lower()  # type: ignore
        if state not in get_args(FinishType):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return

        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"

        console.info("SwanLab Run is finishing, waiting for logs to flush...")

        # 发送强类型的退出信号
        self._queue.put(FinishEvent())

        # 阻塞主线程，等待后台队列消费完毕
        if self._background_logger.is_alive():
            self._background_logger.join()

        console.info(f"Run finished with state: {state}")

    def _define_scalar(
        self,
        key: str,
        name: Optional[str],
        color: Optional[str],
        x_axis: ScalarXAxisType = "_step",
        system: bool = False,
        chart_name: Optional[str] = None,
        chart_index: Optional[str] = None,
    ):
        """
        定义一个标量指标
        :param key: 指标键
        :param name: 指标名称，默认为列名称
        :param color: 指标颜色
        :param x_axis: x轴，可以是其他的标量，也可以是系统值"_step"或"_relative_time"
        :param system: 是否为系统指标
        :param chart_name: 图表名称，默认为列名称
        :param chart_index: 图表索引
        """

        self._queue.put(
            DefineEvent(
                key=key,
                name=name,
                color=color,
                system=system,
                x_axis=x_axis,
                chart_name=chart_name,
                chart=chart_index,
            ),
            block=True,
        )

    # ----------------------------------
    # 消费者核心：统一处理时序与竞态
    # ----------------------------------

    def _background_log(self, flush_timeout: float = 0.5, batch_size: int = 100) -> None:
        """后台单线程
        负责接收并解析事件和罗盘
        """
        batch_records: FlushPayload = []
        # 记录已创建的列，完全避免多线程锁竞争
        _emitted_columns = set()

        while True:
            try:
                # 带超时的阻塞获取，完美平衡吞吐量和延迟
                event = self._queue.get(timeout=flush_timeout)

                # 1. 退出信号
                if isinstance(event, FinishEvent):
                    self._flush(batch_records)
                    break

                # 2. 显式创建列 (Explicit Define)
                elif isinstance(event, DefineEvent):
                    if event.key not in _emitted_columns:
                        col_record = self._dispatch_define(event, key=event.key)
                        batch_records.append(col_record)
                        _emitted_columns.add(event.key)
                    else:
                        console.warning(f"Column '{event.key}' has already been defined, cannot redefine.")

                # 3. 记录数据 (可能触发隐式创建 Implicit Define)
                elif isinstance(event, LogEvent):
                    for key, value in event.data.items():
                        try:
                            # a. 正常解析日志并生成 LogRecord
                            log_record, _ = self._dispatch_parse(
                                value, key=key, timestamp=event.timestamp, step=event.step
                            )
                            # b. 隐式创建检查：如果遇到新 key，自动推断并建列
                            if key not in _emitted_columns:
                                col_record = self._dispatch_define(log_record, key=key)
                                batch_records.append(col_record)
                                _emitted_columns.add(key)
                            # c. 更新上下文中的指标状态，目前仅更新标量指标状态
                            if isinstance(log_record.scalar, ScalarValue):
                                self._metrics.update_scalar(key, log_record.scalar.number)
                            batch_records.append(log_record)
                        except Exception as e:
                            console.error(f"Error when parsing metric '{key}': {e}")

                # 4. 微批处理落盘检查
                if len(batch_records) >= batch_size:
                    self._flush(batch_records)
                    batch_records.clear()

            except queue.Empty:
                # 超时强制刷盘
                if batch_records:
                    self._flush(batch_records)
                    batch_records.clear()
            except Exception as e:
                # 终极防线
                console.error(f"SwanLab background logger thread error: {e}")

    def _flush(self, records: FlushPayload):
        """
        刷盘与回调触发

        1. 将日志和列定义批量写入本地文件
        2. 分发给各个回调

        :param records: 待写入的记录列表
        """
        if not records:
            return
        try:
            # TODO: 批量追加到本地文件
            # 落盘成功后，通过 CallbackManager 异步分发给各个端
            self._callbacker.on_batch_log(records)

        except Exception as e:
            console.error(f"SwanLab failed to write disk or trigger callbacks: {e}")

    # ----------------------------------
    # 多态分发器 (Parsers & Definers)
    # ----------------------------------

    @singledispatchmethod
    def _dispatch_parse(self, value: Any, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """解析标量生成 LogRecord（默认回退：标量）"""
        scalar_value = Scalar.transform(value)
        return LogRecord(key=key, step=step, timestamp=timestamp, scalar=scalar_value), "scalar"

    @_dispatch_parse.register
    def _(self, value: TransformMediaType, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """解析媒体对象生成 LogRecord"""
        cls = value.__class__
        cls_type = cls.type()
        path = self._ctx.media_dir / cls_type
        safe_mkdir(path)
        media_value = cls.transform(key=key, step=step, path=path, content=value)
        record = LogRecord(key=key, step=step, timestamp=timestamp)
        getattr(record, cls_type).CopyFrom(media_value)
        return record, cls_type

    @singledispatchmethod
    def _dispatch_define(self, value: Any, key: str) -> ColumnRecord:
        """
        定义列并生成 ColumnRecord，此时为显式创建列
        """
        # assert isinstance(value, DefineEvent), "Invalid event for column definition."
        ...

    @_dispatch_define.register
    def _(self, value: LogRecord, key: str) -> ColumnRecord:
        """
        解析指标并生成 ColumnRecord
        :param value:
        :param key:
        :return:
        """
        ...
