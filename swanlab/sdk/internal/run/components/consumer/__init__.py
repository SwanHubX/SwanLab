"""
@author: cunyue
@file: consumer.py
@time: 2026/3/13
@description: 后台消费者组件，用于消费运行事件
"""

import math
import queue
import threading
from abc import ABC
from typing import List, Tuple

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord, ScalarValue
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.bus.emitter import RunQueue
from swanlab.sdk.internal.bus.events import (
    ConfigEvent,
    FileSaveEvent,
    LogEvent,
    MetricDefineEvent,
    MetricLogEvent,
)
from swanlab.sdk.internal.context import RunContext, TransformMedia
from swanlab.sdk.internal.pkg import console, helper, safe
from swanlab.sdk.internal.run.components.builder import RecordBuilder
from swanlab.sdk.internal.run.components.resolver import SYSTEM_X_AXES, DefinitionResolver, is_custom_x
from swanlab.sdk.internal.run.transforms import Scalar, echarts

LogBatch = List[LogRecord]
ColumnBatch = List[ColumnRecord]
ScalarBatch = List[ScalarRecord]
MediaBatch = List[MediaRecord]
SaveBatch = List[SaveRecord]
BatchTuple = Tuple[LogBatch, ColumnBatch, ScalarBatch, MediaBatch, SaveBatch]


_EchartsType = (echarts.Base, echarts.Table)


def _is_scalar_value(value: object) -> bool:
    """预判断 value 是否会被 build_scalar_or_media 处理为标量。

    不能调用 build_scalar_or_media：media transform 会往 media_dir 落盘。
    """
    if isinstance(value, TransformMedia):
        return False
    if isinstance(value, list):
        return False
    if isinstance(value, _EchartsType):
        return False
    return True


class ConsumerProtocol(ABC):
    """消费者协议"""

    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        builder: RecordBuilder,
        resolver: DefinitionResolver,
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        _ = (ctx, event_queue, builder, resolver, flush_timeout, batch_size)

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
        resolver: DefinitionResolver,
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        super().__init__(ctx, event_queue, builder, resolver, flush_timeout, batch_size)
        self._ctx = ctx
        self._queue = event_queue
        self._builder = builder
        self._resolver = resolver
        self._core = ctx.core
        self._flush_timeout = flush_timeout
        self._batch_size = batch_size
        self._thread = threading.Thread(target=self._run, name="SwanLab·Specter", daemon=True)
        # 回调器，负责触发回调
        self._callbacker = self._ctx.callbacker

        self._log_batch: LogBatch = []
        self._column_batch: ColumnBatch = []
        self._media_batch: MediaBatch = []
        self._scalar_batch: ScalarBatch = []
        self._save_batch: SaveBatch = []

    def start(self) -> None:
        self._thread.start()

    def join(self) -> None:
        if self._thread.is_alive():
            self._thread.join()

    def stop(self) -> None:
        self._queue.put(_STOP)  # type: ignore[arg-type]

    @property
    def _batch_len(self) -> int:
        return (
            len(self._log_batch)
            + len(self._column_batch)
            + len(self._scalar_batch)
            + len(self._media_batch)
            + len(self._save_batch)
        )

    @property
    def _batch_full(self) -> bool:
        return self._batch_len >= self._batch_size

    @property
    def _batch_empty(self) -> bool:
        return self._batch_len == 0

    def _take_batches(self) -> BatchTuple:
        batches = (
            self._log_batch,
            self._column_batch,
            self._scalar_batch,
            self._media_batch,
            self._save_batch,
        )
        self._log_batch = []
        self._column_batch = []
        self._scalar_batch = []
        self._media_batch = []
        self._save_batch = []
        return batches

    def _restore_batches(
        self,
        log_batch: LogBatch,
        column_batch: ColumnBatch,
        scalar_batch: ScalarBatch,
        media_batch: MediaBatch,
        save_batch: SaveBatch,
    ) -> None:
        # 失败的旧数据插回队头，优先于 flush 期间新进来的数据重试
        if log_batch:
            self._log_batch[:0] = log_batch
        if column_batch:
            self._column_batch[:0] = column_batch
        if scalar_batch:
            self._scalar_batch[:0] = scalar_batch
        if media_batch:
            self._media_batch[:0] = media_batch
        if save_batch:
            self._save_batch[:0] = save_batch

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
        elif isinstance(event, MetricDefineEvent):
            self._resolver.handle_define(event)
        elif isinstance(event, ConfigEvent):
            self._save_batch.append(self._builder.build_config(event))
        elif isinstance(event, LogEvent):
            self._log_batch.append(self._builder.build_log(event))
            console.debug(f"Buffered terminal log record: pending={len(self._log_batch)}", write_to_tty=False)
        elif isinstance(event, FileSaveEvent):
            self._save_batch.append(self._builder.build_save(event))

    def _handle_metric_log(self, event: MetricLogEvent) -> None:
        data = event.data

        # 1. 预扫描：收集显式 scalar 值，更新 cache + 登记 REAL 来源。
        #    仅在登记过 custom X 规则时执行：update_custom_x_cache 以 _custom_x_keys 为门槛，
        #    explicit_scalars 的消费点也全部位于 is_custom_x 分支之后，无 custom X 时本段为死值。
        explicit_scalars: dict[str, float] = {}
        scalar_values: dict[str, ScalarValue] = {}
        if self._resolver.has_custom_x:
            for key, value in data.items():
                if helper.is_system_key(key):
                    continue
                if not _is_scalar_value(value):
                    continue
                with safe.block(message=f"Error when pre-scanning scalar '{key}'", write_to_tty=False):
                    scalar_value = Scalar.transform(value)
                    # transform 结果（含 nan）缓存供第 3 步复用，避免 build 时二次提取 tensor/numpy
                    scalar_values[key] = scalar_value
                    val = scalar_value.number
                    if math.isfinite(val):
                        explicit_scalars[key] = val

            for key, val in explicit_scalars.items():
                self._resolver.update_custom_x_cache(key, val, event.step)

        section_rule_index = self._ctx.config.settings.core.section_rule

        # 2. 计算候选注入（不 commit）：为 custom X 且 event 未含 X 的 Y 缓存最近真实 X 值。
        #    此处只算候选值、不调 try_inject_x——commit 推迟到第 4 步确认有存活 Y 之后，
        #    否则 Y 被去重丢弃后仍会留下孤儿注入点 / 误标 INJECTED 触发误告警（S2）。
        #    仅 step_sync=True（默认）的 Y 会成为注入原因；False 的 Y 不加入候选。
        candidate_x: dict[str, float] = {}
        for key in explicit_scalars:
            concrete = self._resolver.resolve_concrete(key, "SCALAR")
            if not concrete.effective.step_sync:
                continue
            x_axis = concrete.effective.x_axis
            if x_axis in SYSTEM_X_AXES:
                continue
            # 跳过条件看 event.data（用户是否提供过该 key），而非 explicit_scalars：
            # 否则 X 为 nan / 非法值时会被 isfinite 过滤掉，误判为"未提供"而注入，
            # 随后注入值覆盖用户显式值（M1 数据保真 bug）。
            if x_axis in data or x_axis in candidate_x:
                continue
            cached_x = self._resolver.get_custom_x(x_axis)
            if cached_x is None:
                continue
            candidate_x[x_axis] = cached_x

        # 3. 构建 data records + X 值去重，记录“被存活 Y 消费的候选 X”
        consumed_candidate_x: set[str] = set()
        for key, value in data.items():
            with safe.block(message=f"Error when parsing metric '{key}'"):
                cached = scalar_values.get(key)
                if cached is not None:
                    # 预扫描已 transform（同一纯函数、同一输入），直接建 record，
                    # 避免 build_scalar_or_media 对 tensor/numpy 的二次 .item() 提取
                    data_record = Scalar.build_data_record(
                        key=key, step=event.step, timestamp=event.timestamp, data=cached
                    )
                else:
                    data_record, _ = self._builder.build_scalar_or_media(value, key, event.timestamp, event.step)
                if data_record is None:
                    console.warning(f"Metric '{key}' at step {event.step} returned no data, skipped")
                    continue

                if not helper.is_system_key(key):
                    is_scalar = isinstance(data_record, ScalarRecord)
                    metric_class = "SCALAR" if is_scalar else "MEDIA"
                    concrete = self._resolver.resolve_concrete(key, metric_class)
                    col = self._resolver.materialize_column(key, metric_class, data_record.type, section_rule_index)
                    if col is not None:
                        self._column_batch.append(col)

                    # custom X 轴 first-writer-wins：同一 X 值上首次 Y 值为准
                    if is_scalar and is_custom_x(concrete.effective.x_axis):
                        x_axis = concrete.effective.x_axis
                        x_value = explicit_scalars.get(x_axis)
                        used_candidate = False
                        if x_value is None and concrete.effective.step_sync:
                            # True：可用跨 step 的 cache（candidate）。
                            # False：只用本 event 里的显式 X，避免没绑到 X 的 Y 被 cache 误删。
                            x_value = candidate_x.get(x_axis)
                            used_candidate = x_value is not None
                        if x_value is not None and not self._resolver.try_accept_x_value(key, x_value):
                            console.debug(
                                f"Skipping '{key}': duplicate X value {x_value} for '{x_axis}'",
                                write_to_tty=False,
                            )
                            continue
                        # 只有该 Y 存活且确实消费了候选 X，才在后续 commit 注入
                        if used_candidate:
                            consumed_candidate_x.add(x_axis)

                if isinstance(data_record, ScalarRecord):
                    self._scalar_batch.append(data_record)
                else:
                    self._media_batch.append(data_record)

        # 4. 仅为“至少一个存活 Y 消费”的候选 X commit 注入：标 INJECTED + 构建注入 record。
        #    Y 全被去重丢弃时该 X 不在 consumed，从而不产生孤儿注入点、不误标 INJECTED（不触发误告警）。
        for x_axis in candidate_x:
            if x_axis not in consumed_candidate_x:
                continue
            if not self._resolver.try_inject_x(x_axis, event.step):
                continue
            with safe.block(message=f"Error when injecting custom X '{x_axis}'"):
                injected_record, _ = self._builder.build_scalar_or_media(
                    candidate_x[x_axis], x_axis, event.timestamp, event.step
                )
                # 注入值恒为标量（来自 cache 的 float）；isinstance 同时收窄 pyright 类型
                if isinstance(injected_record, ScalarRecord):
                    self._scalar_batch.append(injected_record)

    def _flush(self) -> None:
        if self._batch_empty:
            return
        log_batch, column_batch, scalar_batch, media_batch, save_batch = self._take_batches()
        # 某一步失败时，只回塞"当前未成功提交"的部分，避免重复写入
        # 提交失败时静默显示在 debug 日志中，不打印到控制台
        with safe.block(
            message="Error when flushing batch",
            write_to_tty=False,
            on_error=lambda _: self._restore_batches(log_batch, column_batch, scalar_batch, media_batch, save_batch),
        ):
            if log_batch:
                console.debug(f"Flushing log records to core: count={len(log_batch)}", write_to_tty=False)
                self._core.upsert_logs(log_batch)
                console.debug(f"Core log call completed: count={len(log_batch)}", write_to_tty=False)
                self._callbacker.on_log_flush(log_batch)
                log_batch = []

            if column_batch:
                # 一次 flush 内按 (key, type) coalesce，保留最后一条
                coalesced: dict[tuple[str, int], ColumnRecord] = {}
                for col in column_batch:
                    coalesced[(col.column_key, col.column_type)] = col
                coalesced_list = list(coalesced.values())
                console.debug(f"Flushing column records to core: count={len(coalesced_list)}", write_to_tty=False)
                self._core.upsert_columns(coalesced_list)
                column_batch = []

            if scalar_batch:
                self._core.upsert_scalars(scalar_batch)
                self._callbacker.on_scalar_flush(scalar_batch)
                scalar_batch = []

            if media_batch:
                self._core.upsert_media(media_batch)
                self._callbacker.on_media_flush(media_batch)
                media_batch = []

            if save_batch:
                self._core.upsert_saves(save_batch)
                save_batch = []
