"""
@author: cunyue
@file: record_builder.py
@time: 2026/3/13
@description: protobuf Record 工厂，覆盖所有 record_type
"""

from functools import singledispatchmethod
from typing import Optional, Tuple

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.data.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType, SectionType
from swanlab.proto.swanlab.data.v1.metric_pb2 import MetricRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, RunState
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.proto.swanlab.system.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.sdk.internal.context import RunContext, TransformMediaType
from swanlab.sdk.internal.pkg.fs import safe_mkdir
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.data import MediaTransferType

from .data.transforms import Scalar
from .events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    DefineEvent,
    MetadataEvent,
    ParseResult,
    RequirementsEvent,
    RunStartEvent,
)


class RecordBuilder:
    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        # 由 BackgroundConsumer 单线程调用，无需锁
        self._num: int = 0

    def _wrap(self, **kwargs) -> Record:
        """统一附加 num(自增) + timestamp，返回 Record envelope"""
        self._num += 1
        ts = Timestamp()
        ts.GetCurrentTime()
        return Record(num=self._num, timestamp=ts, **kwargs)

    # ── 用户数据 ──

    @singledispatchmethod
    def build_log(self, value, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """默认回退：标量"""
        scalar_value = Scalar.transform(value)
        metric = MetricRecord(key=key, step=step, timestamp=timestamp, scalar=scalar_value)
        return self._wrap(metric=metric), "scalar"

    @build_log.register
    def _(self, value: TransformMediaType, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """媒体对象"""
        cls = value.__class__
        cls_type = cls.type()
        path = self._ctx.media_dir / cls_type
        safe_mkdir(path)
        media_value = cls.transform(key=key, step=step, path=path, content=value)
        metric = MetricRecord(key=key, step=step, timestamp=timestamp)
        getattr(metric, cls_type).CopyFrom(media_value)
        return self._wrap(metric=metric), cls_type

    def build_column_from_log(self, metric_record: MetricRecord, key: str) -> Record:
        """隐式创建列：从 MetricRecord 推断 ColumnType，并同步 RunMetrics"""
        col_type, section_type = self._infer_column_type(metric_record)
        metrics = self._ctx.metrics
        if col_type == ColumnType.COLUMN_TYPE_FLOAT:
            metrics.define_scalar(key)
        else:
            media_type = self._metric_to_media_type(metric_record)
            metrics.define_media(key, media_type, self._ctx.media_dir / media_type)
        col = ColumnRecord(
            key=key,
            type=col_type,
            section_type=section_type,
            class_=ColumnClass.COLUMN_CLASS_CUSTOM,
        )
        return self._wrap(column=col)

    def build_column_from_define(self, event: DefineEvent) -> Record:
        """显式创建列（DefineEvent），仅支持 scalar"""
        metrics = self._ctx.metrics
        metrics.define_scalar(
            key=event.key,
            name=event.name,
            color=event.color,
            x_axis=event.x_axis,
            system=event.system,
            chart=event.chart,
            chart_name=event.chart_name,
        )
        section_type = SectionType.SECTION_TYPE_SYSTEM if event.system else SectionType.SECTION_TYPE_PUBLIC
        col = ColumnRecord(
            key=event.key,
            type=ColumnType.COLUMN_TYPE_FLOAT,
            class_=ColumnClass.COLUMN_CLASS_CUSTOM,
            section_type=section_type,
            chart_index=event.chart or "",
            chart_name=event.chart_name or "",
            metric_name=event.name or "",
            metric_colors=[event.color] if event.color else [],
        )
        return self._wrap(column=col)

    # ── Run 生命周期 ──

    def build_run(self, event: RunStartEvent) -> Record:
        """构建 RunRecord envelope"""
        return self._wrap(run=event.run_record)

    def build_finish(self, state: FinishType, error: Optional[str] = None) -> Record:
        """构建 FinishRecord envelope"""
        state_map = {
            "success": RunState.RUN_STATE_FINISHED,
            "crashed": RunState.RUN_STATE_CRASHED,
            "aborted": RunState.RUN_STATE_STOPPED,
        }
        ts = Timestamp()
        ts.GetCurrentTime()
        finish = FinishRecord(
            state=state_map.get(state, RunState.RUN_STATE_FINISHED),
            error=error or "",
            finished_at=ts,
        )
        return self._wrap(finish=finish)

    # ── 系统元数据 ──

    def build_config(self, event: ConfigEvent) -> Record:
        """构建 ConfigRecord envelope"""
        return self._wrap(config=event.config_record)

    def build_console(self, event: ConsoleEvent) -> Record:
        """构建 ConsoleRecord envelope"""
        ts = Timestamp()
        ts.GetCurrentTime()
        console_record = ConsoleRecord(line=event.line, stream=event.stream, timestamp=ts)
        return self._wrap(console=console_record)

    def build_metadata(self, event: MetadataEvent) -> Record:
        """构建 MetadataRecord envelope"""
        return self._wrap(metadata=MetadataRecord(path=event.path))

    def build_requirements(self, event: RequirementsEvent) -> Record:
        """构建 RequirementsRecord envelope"""
        return self._wrap(requirements=RequirementsRecord(path=event.path))

    def build_conda(self, event: CondaEvent) -> Record:
        """构建 CondaRecord envelope"""
        return self._wrap(conda=CondaRecord(path=event.path))

    # ── 内部工具 ──

    def _infer_column_type(self, metric: MetricRecord) -> Tuple["ColumnType", "SectionType"]:
        """根据 MetricRecord.value oneof 推断 ColumnType"""
        field_name = metric.WhichOneof("value")
        section_type = SectionType.SECTION_TYPE_PUBLIC
        mapping = {
            "scalar": ColumnType.COLUMN_TYPE_FLOAT,
            "images": ColumnType.COLUMN_TYPE_IMAGE,
            "audios": ColumnType.COLUMN_TYPE_AUDIO,
            "texts": ColumnType.COLUMN_TYPE_TEXT,
            "videos": ColumnType.COLUMN_TYPE_ANY,
            "echarts": ColumnType.COLUMN_TYPE_ANY,
        }
        col_type = mapping.get(field_name, ColumnType.COLUMN_TYPE_UNSPECIFIED)
        return col_type, section_type

    def _metric_to_media_type(self, metric: MetricRecord) -> MediaTransferType:
        """将 MetricRecord oneof 字段名映射为 MediaTransferType"""
        field_name = metric.WhichOneof("value")
        mapping: dict = {
            "images": "image",
            "audios": "audio",
            "texts": "text",
            "videos": "video",
            "echarts": "echarts",
        }
        return mapping.get(field_name, "image")
