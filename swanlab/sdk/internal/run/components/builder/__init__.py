"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: protobuf Record 工厂组件，覆盖除了生命周期以外的所有 record_type
"""

from functools import singledispatchmethod
from typing import List, Type

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import (
    ColumnClass,
    ColumnRecord,
    ColumnType,
    MetricColors,
    SectionType,
)
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.sdk.internal.bus.events import ConfigEvent, ConsoleEvent, ParseResult, ScalarDefineEvent
from swanlab.sdk.internal.context import RunContext, TransformMedia
from swanlab.sdk.internal.context.transformer import TransformData
from swanlab.sdk.internal.pkg import adapter, fs
from swanlab.sdk.internal.run.transforms import Scalar


class RecordBuilder:
    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        # 由 BackgroundConsumer 单线程调用，无需锁
        self._num: int = 0

    # ── 用户数据 ──

    @singledispatchmethod
    def build_log(self, value, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """默认回退：标量"""
        scalar_value = Scalar.transform(value)
        return Scalar.build_data_record(key=key, step=step, timestamp=timestamp, data=scalar_value), Scalar

    @build_log.register(list)
    def _(self, value: List[TransformMedia], key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """媒体对象数组
        dispatch 并不能识别每个数组元素的类型，因此还需手动检查
        """
        if not value or not isinstance(value[0], TransformMedia):
            raise TypeError("List must contain TransformMediaType objects")
        cls = value[0].__class__
        if not all(isinstance(item, cls) for item in value):
            raise TypeError(f"All items in the list must be of the same type {cls.__name__}, got mixed types.")
        path = self._ctx.media_dir / adapter.medium[cls.column_type()]
        fs.safe_mkdir(path)
        values = [item.transform(step=step, path=path) for item in value]
        return cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values), cls

    @build_log.register(TransformMedia)
    def _(self, value: TransformMedia, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """将单个 TransformMediaType 转换为 MediaRecord"""
        cls = value.__class__
        path = self._ctx.media_dir / adapter.medium[cls.column_type()]
        fs.safe_mkdir(path)
        values = [value.transform(step=step, path=path)]
        return cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values), cls

    def build_column_from_log(self, cls: Type[TransformData], key: str) -> ColumnRecord:
        """隐式创建列：从 TransformType 推断 ColumnType，并同步 RunMetrics"""
        column_record = ColumnRecord(
            column_key=key,
            column_type=cls.column_type(),
            column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
            # 自动创建的section默认为公共section
            section_type=SectionType.SECTION_TYPE_PUBLIC,
        )
        col_type = column_record.column_type
        metrics = self._ctx.metrics
        if issubclass(cls, TransformMedia):
            media_type_str = adapter.medium[col_type]
            metrics.define_media(key, column_record, self._ctx.media_dir / media_type_str)
        else:
            metrics.define_scalar(key=key, column=column_record)

        return column_record

    def build_column_from_scalar_define(self, event: ScalarDefineEvent) -> ColumnRecord:
        """显式创建标量列（DefineEvent）"""
        metrics = self._ctx.metrics
        section_type = SectionType.SECTION_TYPE_SYSTEM if event.system else SectionType.SECTION_TYPE_PUBLIC
        col = ColumnRecord(
            column_key=event.key,
            column_type=ColumnType.COLUMN_TYPE_SCALAR,
            column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
            section_name=event.chart_name or "",
            section_type=section_type,
            chart_index=event.chart or "",
            chart_name=event.chart_name or "",
            metric_name=event.name or "",
            metric_colors=MetricColors(light=event.color, dark=event.color) if event.color else None,
        )
        metrics.define_scalar(key=event.key, column=col)
        return col

    # ── 系统元数据 ──
    @staticmethod
    def build_config(event: ConfigEvent) -> ConfigRecord:
        """构建 ConfigRecord envelope"""
        return ConfigRecord(update_type=event.update, timestamp=event.timestamp)

    @staticmethod
    def build_console(event: ConsoleEvent) -> ConsoleRecord:
        """构建 ConsoleRecord envelope"""
        return ConsoleRecord(line=event.line, stream=event.stream, timestamp=event.timestamp)
