"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: protobuf Record 工厂组件，覆盖除了生命周期以外的所有 record_type
"""

from functools import singledispatchmethod
from typing import Optional, Type

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import (
    ColumnClass,
    ColumnRecord,
    ColumnType,
    MetricColors,
    SectionType,
)
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.sdk.internal.bus.events import ConfigEvent, ConsoleEvent, ParseResult, ScalarDefineEvent
from swanlab.sdk.internal.context import RunContext, TransformData, TransformMedia
from swanlab.sdk.internal.pkg import adapter, console, fs
from swanlab.sdk.internal.run.transforms import ECharts, Scalar, echarts

_EchartsType = (echarts.Base, echarts.Table)


class RecordBuilder:
    _MEDIA_MAX_SIZE = 10 * 1024**2  # 10 MB
    _MEDIA_MAX_LENGTH = 108

    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        # 由 BackgroundConsumer 单线程调用，无需锁
        self._num: int = 0

    def _ensure_media_size(self, record: MediaRecord) -> Optional[MediaRecord]:
        """
        确保媒体记录长度不超过限制，否则截断
        确保媒体记录每一项大小不超过限制
        """
        items = list(record.value.items)
        # 1. 过滤掉超过单条大小限制的项
        valid_items = []
        for item in items:
            if item.size > self._MEDIA_MAX_SIZE:
                console.warning(
                    f"Media '{item.filename}' in key '{record.key}' step {record.step} "
                    f"exceeds size limit ({item.size} > {self._MEDIA_MAX_SIZE} bytes), dropped"
                )
                continue
            valid_items.append(item)
        # 2. 全部被过滤则返回 None
        if not valid_items:
            return None
        # 3. 超过长度限制则截断
        if len(valid_items) > self._MEDIA_MAX_LENGTH:
            console.warning(
                f"Media record key '{record.key}' step {record.step} "
                f"has {len(valid_items)} items, exceeding limit {self._MEDIA_MAX_LENGTH}, truncated"
            )
            valid_items = valid_items[: self._MEDIA_MAX_LENGTH]
        # 如果没有变化，直接返回原 record
        if len(valid_items) == len(items) and all(a is b for a, b in zip(valid_items, items)):
            return record
        new_record = MediaRecord()
        new_record.CopyFrom(record)
        del new_record.value.items[:]
        new_record.value.items.extend(valid_items)
        return new_record

    # ── 用户数据 ──

    @singledispatchmethod
    def build_log(self, value, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """默认回退：swanlab.echarts 图表自动包装为 ECharts，否则按标量处理"""
        # 如果是 swanlab.echarts 图表，则自动包装为 ECharts，按照媒体处理
        if isinstance(value, _EchartsType):
            return self.build_log(ECharts(value), key, timestamp, step)
        # 否则按标量处理
        scalar_value = Scalar.transform(value)
        return Scalar.build_data_record(key=key, step=step, timestamp=timestamp, data=scalar_value), Scalar

    @build_log.register(list)
    def _(self, value: list, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """媒体对象数组
        dispatch 并不能识别每个数组元素的类型，因此还需手动检查；
        swanlab.echarts 图表对象会被自动包装为 ECharts
        """
        if not value:
            raise TypeError("List must not be empty")
        # 1. pyecharts 图表对象自动包装为 ECharts
        if not isinstance(value[0], TransformMedia):
            if isinstance(value[0], _EchartsType):
                value = [ECharts(item) for item in value]
            else:
                raise TypeError("List must contain TransformMedia objects or swanlab.echarts chart objects")
        # 2. 确保列表内所有元素类型一致
        cls = value[0].__class__
        if not all(isinstance(item, cls) for item in value):
            raise TypeError(f"All items in the list must be of the same type {cls.__name__}, got mixed types.")
        # 3. 构建媒体记录
        path = self._ctx.media_dir / adapter.medium[cls.column_type()]
        fs.safe_mkdir(path)
        items = [item.transform(step=step, path=path) for item in value]
        media_record = self._ensure_media_size(
            cls.build_data_record(key=key, step=step, timestamp=timestamp, data=items)
        )
        return media_record, cls

    @build_log.register(TransformMedia)
    def _(self, value: TransformMedia, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """将单个 TransformMediaType 转换为 MediaRecord"""
        cls = value.__class__
        path = self._ctx.media_dir / adapter.medium[cls.column_type()]
        fs.safe_mkdir(path)
        values = [value.transform(step=step, path=path)]
        media_record = self._ensure_media_size(
            cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values)
        )
        return media_record, cls

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
