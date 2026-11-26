"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13
@description: protobuf Record 工厂组件，覆盖除了生命周期以外的所有 record_type
"""

from functools import singledispatchmethod
from pathlib import Path
from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarValue, ScoreRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.bus.events import ConfigEvent, FileSaveEvent, LogEvent, ParseResult, ScoreLogEvent
from swanlab.sdk.internal.context import RunContext, TransformMedia
from swanlab.sdk.internal.core_python.metrics.score import ScoreClock
from swanlab.sdk.internal.core_python.store import LEVELDBLOG_HEADER_LEN, DataStoreError, DataStoreReader
from swanlab.sdk.internal.pkg import adapter, console, fs
from swanlab.sdk.internal.run.components.config.writer import dump_config
from swanlab.sdk.internal.run.transforms import ECharts, Histogram, Scalar, echarts

_EchartsType = (echarts.Base, echarts.Table)

# 会被 build_scalar_or_media 处理为媒体（非标量）的类型
_NON_SCALAR_TYPES = (list, TransformMedia, Histogram, *_EchartsType)


def is_scalar_value(value: object) -> bool:
    """判断是否会被 build_scalar_or_media 处理为标量

    在不触发 media transform 落盘的前提下决定是否预 transform。
    判定必须与 build_scalar_or_media 的 dispatch 结果一致：显式注册类型（TransformMedia
    各媒体类与 Histogram 直方图）走非标量分支，默认分支中的 echarts 类型会被包装为
    ECharts 媒体，其余按标量处理。
    """
    return not isinstance(value, _NON_SCALAR_TYPES)


class RecordBuilder:
    _MEDIA_MAX_LENGTH = 108

    _SCORE_CLOCK_STATE = "score_clock.json"

    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        # 分数事件时间分配器：按 (evalId, dataset, key) 槽位严格递增，水位落盘以覆盖重启。
        # 懒创建：从不写分数的 run 不会产生多余状态文件。
        self._score_clock_instance: Optional[ScoreClock] = None
        # 由 BackgroundConsumer 单线程调用，无需锁
        self._num: int = 0
        # 所有媒体类型共用同一单文件上限，由 settings 统一配置
        self._media_max_size: int = ctx.config.settings.core.media_max_size

    def _resolve_media_dir(self, column_type: ColumnType) -> Optional[Path]:
        """解析媒体存储目录。

        skip_store 下返回 None：不创建 media 目录，由 transform 将内容写入 MediaItem.payload。
        """
        if self._ctx.config.settings.core.skip_store:
            return None
        path = self._ctx.media_dir / adapter.medium[column_type]
        fs.safe_mkdir(path)
        return path

    @property
    def _score_clock(self) -> ScoreClock:
        if self._score_clock_instance is None:
            run_dir = getattr(self._ctx, "run_dir", None)
            state_path = Path(run_dir) / self._SCORE_CLOCK_STATE if run_dir is not None else None
            self._score_clock_instance = ScoreClock(state_path, recover=self._recover_score_records)
        return self._score_clock_instance

    def _recover_score_records(self) -> list[ScoreRecord]:
        """从本地数据存储读出已写入的分数记录，用于重建版本水位。

        记录先于上传落盘（core 侧先写 store 再 transport），所以它是进程崩溃或账本
        丢失后唯一可信的事实来源。
        """
        run_file = getattr(self._ctx, "run_file", None)
        if run_file is None or not Path(run_file).exists():
            return []
        size = Path(run_file).stat().st_size
        reader = DataStoreReader()
        try:
            reader.open(str(run_file))
        except (DataStoreError, AssertionError):
            # 文件里连头部都还没落盘：本进程刚创建它，没有任何历史记录。文件不是空的
            # 却读不出头部，才是真的损坏，交给调用侧按"不安全"处理。
            if size < LEVELDBLOG_HEADER_LEN:
                return []
            raise
        records: list[ScoreRecord] = []
        try:
            for payload in reader:
                record = Record()
                record.ParseFromString(payload)
                if record.HasField("score"):
                    records.append(record.score)
        except (DataStoreError, AssertionError):
            # 尾部可能是本进程正在写的记录：已读到的部分仍然可信，内存水位也覆盖了它
            pass
        return records

    def _ensure_media_size(self, record: MediaRecord) -> Optional[MediaRecord]:
        """
        确保媒体记录长度不超过限制，否则截断
        确保媒体记录每一项大小不超过限制
        """
        items = list(record.value.items)
        # 1. 过滤掉超过单条大小限制的项
        valid_items = []
        limit = self._media_max_size
        for item in items:
            if item.size > limit:
                console.warning(
                    f"Media '{item.filename}' in key '{record.key}' step {record.step} "
                    f"exceeds size limit ({item.size} > {limit} bytes), dropped"
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
    def build_scalar_or_media(self, value, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """默认回退：swanlab.echarts 图表自动包装为 ECharts，否则按标量处理"""
        # 如果是 swanlab.echarts 图表，则自动包装为 ECharts，按照媒体处理
        if isinstance(value, _EchartsType):
            return self.build_scalar_or_media(ECharts(value), key, timestamp, step)
        # 否则按标量处理
        scalar_value = Scalar.transform(value)
        return Scalar.build_data_record(key=key, step=step, timestamp=timestamp, data=scalar_value), Scalar

    @build_scalar_or_media.register(list)
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
        path = self._resolve_media_dir(cls.column_type())
        items = [item.transform(step=step, path=path) for item in value]
        media_record = self._ensure_media_size(
            cls.build_data_record(key=key, step=step, timestamp=timestamp, data=items)
        )
        return media_record, cls

    @build_scalar_or_media.register(TransformMedia)
    def _(self, value: TransformMedia, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """将单个 TransformMediaType 转换为 MediaRecord"""
        cls = value.__class__
        path = self._resolve_media_dir(cls.column_type())
        values = [value.transform(step=step, path=path)]
        media_record = self._ensure_media_size(
            cls.build_data_record(key=key, step=step, timestamp=timestamp, data=values)
        )
        return media_record, cls

    @build_scalar_or_media.register(Histogram)
    def _(self, value: Histogram, key: str, timestamp: Timestamp, step: int) -> ParseResult:
        """将 Histogram 对象转换为 HistogramRecord"""
        # 类型/形状校验失败已在构造时记录日志，这里直接跳过该 key
        if not value.is_valid:
            return None, Histogram
        data_record = Histogram.build_data_record(key=key, step=step, timestamp=timestamp, data=value.transform())
        return data_record, Histogram

    # ── 系统元数据 ──
    def build_score(self, event: ScoreLogEvent) -> list[ScoreRecord]:
        """把一次 log_score 调用展开为若干条 ScoreRecord。

        值域已在 API 层校验（有限浮点数或 NaN），这里不再做归一化，避免二次转换把
        同一个输入解释成另一个值。时间戳按槽位分配一次（ScoreClock），此后随 record
        持久化，重放不再改动。
        """
        records: list[ScoreRecord] = []
        for key, number in event.scores.items():
            record = ScoreRecord(
                eval_id=event.eval_id,
                dataset=event.dataset,
                key=key,
                value=ScalarValue(number=number),
            )
            record.timestamp.CopyFrom(self._score_clock.next_timestamp((event.eval_id, event.dataset, key)))
            records.append(record)
        return records

    def build_config(self, event: ConfigEvent) -> SaveRecord:
        """构建 Config Save envelope。

        skip_store 下 config 不落盘，内容按落盘同款 YAML 编码填入 payload；
        未启用时 payload 为 None，Core 按 source_path 读取 config.yaml。
        """
        if self._ctx.config.settings.core.skip_store:
            return SaveRecord(
                name="config",
                type=SaveType.SAVE_TYPE_CONFIG,
                payload=dump_config(event.content).encode("utf-8"),
            )
        return SaveRecord(name="config", source_path=event.path.absolute().as_posix(), type=SaveType.SAVE_TYPE_CONFIG)

    @staticmethod
    def build_log(event: LogEvent) -> LogRecord:
        """构建 LogRecord envelope"""
        return LogRecord(line=event.line, level=event.level, timestamp=event.timestamp)

    # ── 文件保存 (Save) 数据 ──
    @classmethod
    def build_save(cls, event: FileSaveEvent) -> SaveRecord:
        """将 FileSaveEvent 转为 proto SaveRecord"""
        return SaveRecord(name=event.name, source_path=event.source_path, policy=adapter.policy[event.policy])
