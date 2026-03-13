"""
@author: cunyue
@file: consumer.py
@time: 2026/3/13
@description: SwanLab 后台消费者线程
"""

import queue
import threading
from typing import Set

from swanlab.sdk.internal.bus.emitter import RunQueue
from swanlab.sdk.internal.bus.events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    FlushPayload,
    MetadataEvent,
    MetricDefineEvent,
    MetricLogEvent,
    RequirementsEvent,
    RunFinishEvent,
    RunStartEvent,
)
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core import CoreProtocol
from swanlab.sdk.internal.pkg import console

from .record_builder import RecordBuilder


class BackgroundConsumer:
    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        builder: RecordBuilder,
        core: CoreProtocol,
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        self._ctx = ctx
        self._queue = event_queue
        self._builder = builder
        self._core = core
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

    def _run(self) -> None:
        mode = self._ctx.config.settings.mode
        self._core.startup(cloud=mode == "cloud", persistence=mode != "disabled")
        batch: FlushPayload = []
        # 记录已创建的列，完全避免多线程锁竞争
        _emitted_columns: Set[str] = set()

        while True:
            try:
                # 带超时的阻塞获取，平衡吞吐量和延迟
                event = self._queue.get(timeout=self._flush_timeout)

                # 1. 退出信号
                if isinstance(event, RunFinishEvent):
                    batch.append(self._builder.build_finish(event))
                    self._flush(batch)
                    self._core.shutdown()
                    break

                # 2. 记录数据（可能触发隐式创建 Implicit Define）
                elif isinstance(event, MetricLogEvent):
                    for key, value in event.data.items():
                        try:
                            record, data_type = self._builder.build_log(value, key, event.timestamp, event.step)
                            if key not in _emitted_columns:
                                batch.append(self._builder.build_column_from_log(record.metric, key))
                                _emitted_columns.add(key)
                            if data_type == "scalar":
                                self._metrics.update_scalar(key, record.metric.scalar.number)
                            batch.append(record)
                        except Exception as e:
                            console.error(f"Error when parsing metric '{key}': {e}")

                # 3. 显式创建列 (Explicit Define)
                elif isinstance(event, MetricDefineEvent):
                    if event.key not in _emitted_columns:
                        batch.append(self._builder.build_column_from_define(event))
                        _emitted_columns.add(event.key)
                    else:
                        console.warning(f"Column '{event.key}' has already been defined, cannot redefine.")

                # 4. 系统事件
                elif isinstance(event, RunStartEvent):
                    batch.append(self._builder.build_run(event))
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

                # 5. 微批处理落盘检查
                if len(batch) >= self._batch_size:
                    self._flush(batch)
                    batch.clear()

            except queue.Empty:
                # 超时强制刷盘
                if batch:
                    self._flush(batch)
                    batch.clear()
            except Exception as e:
                # 终极防线
                console.error(f"SwanLab background logger thread error: {e}")

    def _flush(self, records: FlushPayload) -> None:
        """处理一批 Record：持久化、回调、上传等"""
        if not records:
            return
        try:
            self._core.handle_records(records)
            # TODO: 分类数据，根据语义依次触发回调
        except Exception as e:
            console.error(f"SwanLab failed to handle records: {e}")
