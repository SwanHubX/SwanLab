"""
@author: cunyue
@file: record_builder.py
@time: 2026/3/13
@description: protobuf Record 工厂，覆盖所有 record_type
"""

from functools import singledispatchmethod
from typing import List, Type

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType, SectionType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, RunRecord
from swanlab.proto.swanlab.system.v1.console_pb2 import ConsoleRecord
from swanlab.proto.swanlab.system.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.sdk.internal import adapter
from swanlab.sdk.internal.bus.events import (
    CondaEvent,
    ConfigEvent,
    ConsoleEvent,
    MetadataEvent,
    ParseResult,
    RequirementsEvent,
    RunFinishEvent,
    RunStartEvent,
    ScalarDefineEvent,
)
from swanlab.sdk.internal.context import RunContext, TransformMedia
from swanlab.sdk.internal.context.transformer import TransformData
from swanlab.sdk.internal.pkg.fs import safe_mkdir
from swanlab.sdk.internal.run.transforms import Scalar


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
        return self._wrap(
            metric=Scalar.build_data_record(key=key, step=step, timestamp=timestamp, data=scalar_value)
        ), Scalar

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
        path = self._ctx.media_dir / adapter.column_type[cls.column_type()]
        safe_mkdir(path)
        values = [item.transform(step=step, path=path) for item in value]
        return self._wrap(metric=cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values)), cls

    @build_log.register(TransformMedia)
    def _(self, value: TransformMedia, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """将单个 TransformMediaType 转换为 DataRecord"""
        cls = value.__class__
        path = self._ctx.media_dir / adapter.column_type[cls.column_type()]
        safe_mkdir(path)
        values = [value.transform(step=step, path=path)]
        return self._wrap(metric=cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values)), cls

    def build_column_from_log(self, cls: Type[TransformData], key: str) -> Record:
        """隐式创建列：从 TransformType 推断 ColumnType，并同步 RunMetrics"""
        col_type = cls.column_type()
        metrics = self._ctx.metrics
        if issubclass(cls, TransformMedia):
            media_type_str = adapter.column_type[col_type]
            metrics.define_media(key, col_type, self._ctx.media_dir / media_type_str)
        else:
            metrics.define_scalar(key)
        col = ColumnRecord(
            column_key=key,
            column_type=col_type,
            column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
            # 自动创建的section默认为公共section
            section_type=SectionType.SECTION_TYPE_PUBLIC,
        )
        return self._wrap(column=col)

    def build_column_from_scalar_define(self, event: ScalarDefineEvent) -> Record:
        """显式创建标量列（DefineEvent）"""
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
            column_key=event.key,
            column_type=ColumnType.COLUMN_TYPE_FLOAT,
            column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
            section_name=event.chart_name or "",
            section_type=section_type,
            chart_index=event.chart or "",
            chart_name=event.chart_name or "",
            metric_name=event.name or "",
            metric_colors=[event.color, event.color] if event.color else [],
        )
        return self._wrap(column=col)

    # ── Run 生命周期 ──

    def build_run(self, event: RunStartEvent) -> Record:
        """构建 RunRecord envelope"""
        settings = self._ctx.config.settings
        run_record = RunRecord(
            project=settings.project.name,
            workspace=settings.project.workspace,
            name=settings.experiment.name,
            color=settings.experiment.color,
            description=settings.experiment.description,
            job_type=settings.experiment.job_type,
            group=settings.experiment.group,
            tags=settings.experiment.tags,
            id=settings.run.id,
            resume=adapter.resume.get(settings.run.resume),
            started_at=event.timestamp,
        )
        return self._wrap(run=run_record)

    def build_finish(self, event: RunFinishEvent) -> Record:
        """构建 FinishRecord envelope"""
        ts = Timestamp()
        ts.GetCurrentTime()
        finish = FinishRecord(
            state=adapter.state.get(event.state),
            error=event.error or "",
            finished_at=ts,
        )
        return self._wrap(finish=finish)

    # ── 系统元数据 ──

    def build_config(self, event: ConfigEvent) -> Record:
        """构建 ConfigRecord envelope"""
        config_record = ConfigRecord(update_type=event.update, timestamp=event.timestamp)
        return self._wrap(config=config_record)

    def build_console(self, event: ConsoleEvent) -> Record:
        """构建 ConsoleRecord envelope"""
        console_record = ConsoleRecord(line=event.line, stream=event.stream, timestamp=event.timestamp)
        return self._wrap(console=console_record)

    def build_metadata(self, event: MetadataEvent) -> Record:
        """构建 MetadataRecord envelope"""
        return self._wrap(metadata=MetadataRecord(timestamp=event.timestamp))

    def build_requirements(self, event: RequirementsEvent) -> Record:
        """构建 RequirementsRecord envelope"""
        return self._wrap(requirements=RequirementsRecord(timestamp=event.timestamp))

    def build_conda(self, event: CondaEvent) -> Record:
        """构建 CondaRecord envelope"""
        return self._wrap(conda=CondaRecord(timestamp=event.timestamp))
