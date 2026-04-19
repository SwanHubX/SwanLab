"""
@author: cunyue
@file: consumer.py
@time: 2026/3/13
@description: SwanLab 后台消费者线程
"""

import queue
import threading
from abc import ABC
from typing import Set

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.sdk.internal.bus.emitter import RunQueue
from swanlab.sdk.internal.bus.events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    FlushPayload,
    MetadataEvent,
    MetricLogEvent,
    RequirementsEvent,
    ScalarDefineEvent,
)
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, safe

from .builder import RecordBuilder


class ConsumerProtocol(ABC):
    """消费者协议"""

    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        builder: RecordBuilder,
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
        builder: RecordBuilder,
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
        # 指标状态，完全避免多线程锁竞争
        self._metrics = self._ctx.metrics
        # 回调器，负责触发回调
        self._callbacker = self._ctx.callbacker

    def start(self) -> None:
        self._thread.start()

    def join(self) -> None:
        if self._thread.is_alive():
            self._thread.join()

    def stop(self) -> None:
        self._queue.put(_STOP)  # type: ignore[arg-type]

    def _run(self) -> None:
        batch: FlushPayload = []
        # 记录已创建的列，完全避免多线程锁竞争
        _emitted_columns: Set[str] = set()

        while True:
            with safe.block(message="SwanLab background logger thread error"):
                try:
                    # 带超时的阻塞获取，平衡吞吐量和延迟
                    event = self._queue.get(timeout=self._flush_timeout)

                    # 1. 退出信号
                    if event is _STOP:
                        self._flush(batch)
                        break
                    # 3. 记录数据（可能触发隐式创建 Implicit Define）
                    elif isinstance(event, MetricLogEvent):
                        for key, value in event.data.items():
                            with safe.block(message=f"Error when parsing metric '{key}'"):
                                record, cls = self._builder.build_log(value, key, event.timestamp, event.step)
                                if key not in _emitted_columns:
                                    batch.append(self._builder.build_column_from_log(cls, key))
                                    _emitted_columns.add(key)
                                if cls.column_type() == ColumnType.COLUMN_TYPE_FLOAT:
                                    self._metrics.update_scalar(key, record.metric.scalar.number)
                                batch.append(record)
                    # 4. 显式创建列 (Explicit Define)
                    elif isinstance(event, ScalarDefineEvent):
                        if event.key not in _emitted_columns:
                            batch.append(self._builder.build_column_from_scalar_define(event))
                            _emitted_columns.add(event.key)
                        else:
                            console.warning(f"Column '{event.key}' has already been defined, cannot redefine.")
                    # 5. 系统事件
                    elif isinstance(event, ConfigEvent):
                        batch.append(self._builder.build_config(event))
                    elif isinstance(event, ConsoleEvent):
                        batch.append(self._builder.build_console(event))
                    elif isinstance(event, MetadataEvent):
                        batch.append(self._builder.build_metadata(event))
                    elif isinstance(event, RequirementsEvent):
                        batch.append(self._builder.build_requirements(event))
                    elif isinstance(event, CondaEvent):
                        batch.append(self._builder.build_conda(event))

                    # 6. 微批处理落盘检查
                    if len(batch) >= self._batch_size:
                        self._flush(batch)
                        batch.clear()
                except queue.Empty:
                    # 超时强制刷盘
                    if batch:
                        self._flush(batch)
                        batch.clear()

    def _flush(self, records: FlushPayload) -> None:
        """处理一批 Record：持久化、回调、上传等"""
        if not records:
            return
        with safe.block(message="SwanLab failed to handle records"):
            self._core.publish(records)
