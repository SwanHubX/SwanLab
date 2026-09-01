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
from typing import Any, Dict, List, Tuple, Union

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
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, helper, safe
from swanlab.sdk.internal.run.transforms import Scalar

from .builder import RecordBuilder, is_scalar_value
from .resolver import DefinitionResolver

LogBatch = List[LogRecord]
ColumnBatch = List[ColumnRecord]
ScalarBatch = List[ScalarRecord]
MediaBatch = List[MediaRecord]
SaveBatch = List[SaveRecord]
BatchTuple = Tuple[LogBatch, ColumnBatch, ScalarBatch, MediaBatch, SaveBatch]


class ConsumerProtocol(ABC):
    """消费者协议"""

    def __init__(
        self,
        ctx: RunContext,
        event_queue: RunQueue,
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        _ = (ctx, event_queue, flush_timeout, batch_size)

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
        flush_timeout: float = 0.5,
        batch_size: int = 100,
    ):
        super().__init__(ctx, event_queue, flush_timeout, batch_size)
        self._ctx = ctx
        self._queue = event_queue
        # builder / resolver 均由 consumer 独占创建与维护，所有状态访问都发生在本线程内
        self._builder = RecordBuilder(ctx)
        self._resolver = DefinitionResolver()
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
        # 线程退出时清理状态
        self._resolver.clear()

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
        """处理指标日志事件：构建 scalar/media 记录；调用过 define_metric 时额外物化列定义，
        并为 custom X 轴指标做跨 step 的 X/Y 对齐。

        按是否登记过 define rule 分流：
        - 未调用过 define_metric → _log_plain：纯构建，列由 core 收到数据后自动创建；
        - 存在 define rule → _log_define：物化 define 过的列 + custom X 对齐机制。
        """
        # 1. 过滤用户伪造的系统前缀 key
        raw_data = event.data
        data: Dict[str, Any] = {}
        for key, value in raw_data.items():
            if helper.is_system_key(key):
                # 系统列名不允许用户伪造，避免覆盖 swanlab 内部列定义
                # 这种情况比较少见，并且通常为有意为之，因此打印warning而非debug
                console.warning(f"Metric '{key}' at step {event.step} is a system key, skipped")
                continue
            data[key] = value

        if self._resolver.has_rules:
            self._log_define(event, data)
        else:
            self._log_plain(event, data)

    def _log_plain(self, event: MetricLogEvent, data: Dict[str, Any]) -> None:
        """无 define 路径：纯构建 data record，不产出任何 ColumnRecord——列由 core
        收到数据后自动创建（与未引入 define_metric 时的行为一致）。

        仍调用 resolve_concrete 登记 automatic 状态：钉住 "key 首次 log 后 define
        不再生效" 的契约，之后到达的 define 无法认领已 log 过的 key。
        """
        for key, value in data.items():
            with safe.block(message=f"Error when parsing metric '{key}'"):
                data_record, _ = self._builder.build_scalar_or_media(value, key, event.timestamp, event.step)
                if data_record is None:
                    console.warning(f"Metric '{key}' at step {event.step} returned no data, skipped")
                    continue
                metric_class = "SCALAR" if isinstance(data_record, ScalarRecord) else "MEDIA"
                self._resolver.resolve_concrete(key, metric_class)  # 仅登记钉住，不物化
                self._append_data_record(data_record)

    def _log_define(self, event: MetricLogEvent, data: Dict[str, Any]) -> None:
        """有 define 路径：预扫描/候选注入（存在 custom X 规则时）→ 构建 + 物化 define
        过的列 + custom X 的 Y 按 X 值去重 → 提交注入。
        """
        # 2. 预扫描，仅在登记过 custom X 规则时执行：收集显式 scalar 值并更新 custom X 缓存
        explicit_scalars: dict[str, float] = {}
        scalar_values: dict[str, ScalarValue] = {}
        if self._resolver.has_custom_x:
            for key, value in data.items():
                if not is_scalar_value(value):
                    continue
                # 预扫描失败不单独报 error：第 3 步 build 时同一值会再次 transform，
                # 届时以完整上下文上报；此处仅留 debug 级别
                with safe.block(message=f"Error when pre-scanning scalar '{key}'", level="debug", write_to_tty=False):
                    scalar_value = Scalar.transform(value)
                    # transform 结果（含 nan）缓存供第 4 步复用，避免 build 时二次提取 tensor/numpy
                    scalar_values[key] = scalar_value
                    val = scalar_value.number
                    if math.isfinite(val):
                        explicit_scalars[key] = val
            # 收集完毕，更新 custom X 缓存，此时搜集仅显式标量值，避免 media/nan 注入
            for key, val in explicit_scalars.items():
                self._resolver.update_custom_x_cache(key, val, event.step)

        # 3. 计算候选注入，为 custom X 且 event 未含 X 的 Y 生成最近真实 X 值
        # candidate_x 将是最终上报的注入 X 值
        candidate_x: dict[str, float] = {}
        for key in explicit_scalars:
            concrete = self._resolver.resolve_concrete(key, "SCALAR")
            # 3.1 显式 step_sync=False 的 Y 不参与注入
            if not concrete.effective.step_sync:
                continue
            # 3.2 如果x轴为系统内部step，则不参与注入
            x_axis = concrete.effective.x_axis
            if not helper.is_custom_x(x_axis):
                continue
            # 3.3 如果用户在本次 event 里显式 log 了 X 值，则不注入
            if x_axis in data or x_axis in candidate_x:
                continue
            # 3.4 从 cache 取最近真实 X 值；X 尚无任何真实值时注入默认 0，
            # 使 X 序列从首个绑定 Y 的 step 起出现，而非等到 X 自身首次 log
            cached_x = self._resolver.get_custom_x(x_axis)
            if cached_x is None:
                cached_x = 0.0
            candidate_x[x_axis] = cached_x

        # 4. 把本次 log 的每个 key 变成 record，包括：
        # - 首次物化 ColumnRecord
        # - 构建 ScalarRecord / MediaRecord
        # - custom X 标量按 X 值去重，记录被存活 Y 消费的候选 X
        consumed_candidate_x: set[str] = set()
        for key, value in data.items():
            with safe.block(message=f"Error when parsing metric '{key}'"):
                # 4.1 复用预扫描的 transform 结果构建 record，避免对 tensor/numpy 的二次 .item() 提取
                cached = scalar_values.get(key)
                if cached is not None:
                    data_record = Scalar.build_data_record(
                        key=key, step=event.step, timestamp=event.timestamp, data=cached
                    )
                else:
                    data_record, _ = self._builder.build_scalar_or_media(value, key, event.timestamp, event.step)
                if data_record is None:
                    console.warning(f"Metric '{key}' at step {event.step} returned no data, skipped")
                    continue

                # 4.2 物化 ColumnRecord，resolver侧保证仅对define过的 key 物化一次，后续 log 不再重复发 ColumnRecord
                # 对于未define的key，物化行为依旧交给core自动处理
                is_scalar = isinstance(data_record, ScalarRecord)
                metric_class = "SCALAR" if is_scalar else "MEDIA"
                concrete = self._resolver.resolve_concrete(key, metric_class)
                col = self._resolver.materialize_column(key, metric_class, data_record.type)
                if col is not None:
                    self._column_batch.append(col)

                # 4.3 custom X 轴 first-writer-wins：同一 X 值上首次 Y 值为准
                x_axis = concrete.effective.x_axis
                if is_scalar and helper.is_custom_x(x_axis):
                    x_value = explicit_scalars.get(x_axis)
                    used_candidate = False
                    # 进行 step 对齐，如果当前log没有显式 log X 值，则尝试从候选注入中取最近的真实 X 值
                    if x_value is None and concrete.effective.step_sync:
                        x_value = candidate_x.get(x_axis)
                        used_candidate = x_value is not None
                    # 对本次提交的(x, y)进行x轴值的去重，比如(1,2)、(1,3)，仅保留第一条(1,2)，第二条(1,3)会被丢弃
                    if x_value is not None and not self._resolver.try_accept_x_value(key, x_value):
                        console.debug(
                            f"Skipping '{key}': duplicate X value {x_value} for '{x_axis}'",
                            write_to_tty=False,
                            write_to_file=True,
                        )
                        continue
                    if used_candidate:
                        consumed_candidate_x.add(x_axis)

                self._append_data_record(data_record)

        # 5. 提交注入X值：只有当某个 Y 真正用上了缓存 X 值、且该 Y 没被去重丢弃时，才为这个 X 构建 record 并标记本 step 已注入。
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

    def _append_data_record(self, data_record: Union[ScalarRecord, MediaRecord]) -> None:
        """按类型把 data record 追加到对应批次。"""
        if isinstance(data_record, ScalarRecord):
            self._scalar_batch.append(data_record)
        else:
            self._media_batch.append(data_record)

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
