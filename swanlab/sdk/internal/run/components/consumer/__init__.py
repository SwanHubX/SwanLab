"""
@author: cunyue
@file: consumer.py
@time: 2026/3/13
@description: 后台消费者组件，用于消费运行事件
"""

import queue
import threading
from abc import ABC
from typing import TYPE_CHECKING, List

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import DataRecord
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.sdk.internal.bus.emitter import RunQueue
from swanlab.sdk.internal.bus.events import ConfigEvent, ConsoleEvent, MetricLogEvent, ScalarDefineEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.run.components.builder import RecordBuilder

ConfigBatch = List[ConfigRecord]
ConsoleBatch = List[ConsoleRecord]
ColumnBatch = List[ColumnRecord]
DataBatch = List[DataRecord]


class ConsumerProtocol(ABC):
    """消费者协议"""

    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        builder: "RecordBuilder",
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ): ...

    def start(self) -> None: ...

    def stop(self) -> None: ...

    def join(self) -> None: ...


class _StopEvent:
    pass


_STOP = _StopEvent()


class BackgroundConsumer(ConsumerProtocol):
    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        builder: "RecordBuilder",
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        super().__init__(ctx, event_queue, builder, flush_timeout, batch_size)
        self._ctx = ctx
        self._queue = event_queue
        self._builder = builder
        self._core = ctx.core
        self._flush_timeout = flush_timeout
        self._batch_size = batch_size
        self._thread = threading.Thread(target=self._run, name="SwanLab·Specter", daemon=True)
        # 指标状态
        self._metrics = self._ctx.metrics
        # 回调器，负责触发回调
        self._callbacker = self._ctx.callbacker

        self._config_batch: ConfigBatch = []
        self._console_batch: ConsoleBatch = []
        self._column_batch: ColumnBatch = []
        self._data_batch: DataBatch = []

    def start(self) -> None:
        self._thread.start()

    def join(self) -> None:
        if self._thread.is_alive():
            self._thread.join()

    def stop(self) -> None:
        self._queue.put(_STOP)  # type: ignore[arg-type]

    @property
    def _batch_full(self):
        l = len(self._config_batch) + len(self._console_batch) + len(self._column_batch) + len(self._data_batch)  # noqa: E741
        return l >= self._batch_size

    @property
    def _batch_empty(self):
        return not self._config_batch and not self._console_batch and not self._column_batch and not self._data_batch

    def _run(self) -> None:
        while True:
            with safe.block(message="SwanLab background logger thread error"):
                try:
                    # 带超时的阻塞获取，平衡吞吐量和延迟
                    event = self._queue.get(timeout=self._flush_timeout)

                    # 1. 退出信号
                    if event is _STOP:
                        self._flush()
                        break
                    # 3. 记录数据（可能触发隐式创建 Implicit Define）
                    elif isinstance(event, MetricLogEvent):
                        for key, value in event.data.items():
                            with safe.block(message=f"Error when parsing metric '{key}'"):
                                data_record, cls = self._builder.build_log(value, key, event.timestamp, event.step)
                                if not self._metrics.has(key, cls.column_type()):
                                    this_column = self._builder.build_column_from_log(cls, key)
                                    self._column_batch.append(this_column)
                                    self._metrics.define_scalar(key, this_column, value=value)
                                if cls.column_type() == ColumnType.COLUMN_TYPE_FLOAT:
                                    self._metrics.update_scalar(key, data_record.scalar.number)
                                self._data_batch.append(data_record)
                    # 4. 显式创建标量列 (Explicit Define)
                    elif isinstance(event, ScalarDefineEvent):
                        try:
                            if not self._metrics.has(event.key, ColumnType.COLUMN_TYPE_FLOAT):
                                this_column = self._builder.build_column_from_scalar_define(event)
                                self._column_batch.append(this_column)
                                self._metrics.define_scalar(event.key, this_column, value=None)
                        except TypeError:
                            console.warning(f"Scalar Column '{event.key}' has already been defined, cannot redefine.")
                    # 5. 系统事件
                    elif isinstance(event, ConfigEvent):
                        self._config_batch.append(self._builder.build_config(event))
                    elif isinstance(event, ConsoleEvent):
                        self._console_batch.append(self._builder.build_console(event))
                    # 6. 微批处理落盘检查
                    if self._batch_full:
                        self._flush()
                except queue.Empty:
                    # 超时强制刷盘
                    if self._batch_empty:
                        self._flush()

    def _flush(self) -> None:
        """处理一批 Record：持久化、回调、上传等"""
        with safe.block(message="SwanLab failed to flush config batch"):
            if self._config_batch:
                self._core.upsert_configs(self._config_batch)
        with safe.block(message="SwanLab failed to flush console batch"):
            if self._console_batch:
                self._core.upsert_consoles(self._console_batch)
        with safe.block(message="SwanLab failed to flush column batch"):
            if self._column_batch:
                self._core.upsert_columns(self._column_batch)
        with safe.block(message="SwanLab failed to flush data batch"):
            if self._data_batch:
                self._core.upsert_data(self._data_batch)
