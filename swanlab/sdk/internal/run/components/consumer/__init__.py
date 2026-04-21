"""
@author: cunyue
@file: consumer.py
@time: 2026/3/13
@description: 后台消费者组件，用于消费运行事件
"""

import queue
import threading
from abc import ABC
from typing import TYPE_CHECKING, List, Tuple

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
BatchTuple = Tuple[ConfigBatch, ConsoleBatch, ColumnBatch, DataBatch]


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
    def _batch_len(self) -> int:
        return len(self._config_batch) + len(self._console_batch) + len(self._column_batch) + len(self._data_batch)

    @property
    def _batch_full(self) -> bool:
        return self._batch_len >= self._batch_size

    @property
    def _batch_empty(self) -> bool:
        return self._batch_len == 0

    def _take_batches(self) -> BatchTuple:
        batches = (
            self._config_batch,
            self._console_batch,
            self._column_batch,
            self._data_batch,
        )
        self._config_batch = []
        self._console_batch = []
        self._column_batch = []
        self._data_batch = []
        return batches

    def _restore_batches(
        self,
        config_batch: ConfigBatch,
        console_batch: ConsoleBatch,
        column_batch: ColumnBatch,
        data_batch: DataBatch,
    ) -> None:
        # 失败的旧数据插回队头，优先于 flush 期间新进来的数据重试
        if config_batch:
            self._config_batch[:0] = config_batch
        if console_batch:
            self._console_batch[:0] = console_batch
        if column_batch:
            self._column_batch[:0] = column_batch
        if data_batch:
            self._data_batch[:0] = data_batch

    def _run(self) -> None:
        while True:
            with safe.block(message="SwanLab background logger thread error"):
                try:
                    event = self._queue.get(timeout=self._flush_timeout)
                except queue.Empty:
                    # 队列空闲一段时间后，若有积压数据则刷盘
                    if not self._batch_empty:
                        self._flush()
                    continue
                if event is _STOP:
                    self._flush()
                    break
                self._handle_event(event)
                if self._batch_full:
                    self._flush()

    def _handle_event(self, event) -> None:
        if isinstance(event, MetricLogEvent):
            self._handle_metric_log(event)
        elif isinstance(event, ScalarDefineEvent):
            self._handle_scalar_define(event)
        elif isinstance(event, ConfigEvent):
            self._config_batch.append(self._builder.build_config(event))
        elif isinstance(event, ConsoleEvent):
            self._console_batch.append(self._builder.build_console(event))

    def _handle_metric_log(self, event: MetricLogEvent) -> None:
        for key, value in event.data.items():
            with safe.block(message=f"Error when parsing metric '{key}'"):
                data_record, cls = self._builder.build_log(value, key, event.timestamp, event.step)
                defined = self._metrics.ensure_defined_as(key, cls.column_type())
                if not defined:
                    this_column = self._builder.build_column_from_log(cls, key)
                    self._column_batch.append(this_column)
                if cls.column_type() == ColumnType.COLUMN_TYPE_FLOAT:
                    self._metrics.update_scalar(key, data_record.scalar.number)
                self._data_batch.append(data_record)

    def _handle_scalar_define(self, event: ScalarDefineEvent) -> None:
        defined = self._metrics.ensure_defined_as(event.key, ColumnType.COLUMN_TYPE_FLOAT)
        if defined:
            console.warning(f"Scalar Column '{event.key}' has already been defined, skip redefine.")
            return
        this_column = self._builder.build_column_from_scalar_define(event)
        self._column_batch.append(this_column)

    def _flush(self) -> None:
        if self._batch_empty:
            return
        config_batch, console_batch, column_batch, data_batch = self._take_batches()
        # 某一步失败时，只回塞“当前未成功提交”的部分，避免重复写入
        # 提交失败时静默显示在 debug 日志中，不打印到控制台
        with safe.block(
            message="Error when flushing batch",
            write_to_tty=False,
            on_error=lambda _: self._restore_batches(config_batch, console_batch, column_batch, data_batch),
        ):
            if config_batch:
                self._core.upsert_configs(config_batch)
                config_batch = []

            if console_batch:
                self._core.upsert_consoles(console_batch)
                console_batch = []

            if column_batch:
                self._core.upsert_columns(column_batch)
                column_batch = []

            if data_batch:
                self._core.upsert_data(data_batch)
                data_batch = []
