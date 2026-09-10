"""
@author: cunyue
@file: test_builder.py
@time: 2026/4/29
@description: RecordBuilder 单元测试
"""

from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveType
from swanlab.sdk.internal.bus.events import ConfigEvent
from swanlab.sdk.internal.run.components.consumer.builder import _NON_SCALAR_TYPES, RecordBuilder, is_scalar_value
from swanlab.sdk.internal.run.transforms import Text


def _make_media_record(key: str = "test", step: int = 0, items=None) -> MediaRecord:
    """构造一个 MediaRecord，items 为 (filename, size) 元组列表"""
    record = MediaRecord(key=key, step=step)
    if items:
        for filename, size in items:
            record.value.items.append(MediaItem(filename=filename, size=size))
    return record


class TestEnsureMediaSize:
    """RecordBuilder._ensure_media_size 单元测试"""

    # 不依赖 RunContext，直接实例化并覆盖类属性
    @pytest.fixture
    def builder(self):
        b = object.__new__(RecordBuilder)
        b._MEDIA_MAX_SIZE = 100  # 100 bytes  # type: ignore
        b._MEDIA_MAX_LENGTH = 3  # type: ignore
        return b

    def test_all_items_within_limits_returns_original(self, builder):
        """所有项都在限制内，返回原 record"""
        record = _make_media_record(items=[("a.png", 10), ("b.png", 20)])
        result = builder._ensure_media_size(record)
        assert result is record

    def test_single_item_exceeds_size_dropped(self, builder):
        """单条超过大小限制的项被丢弃"""
        record = _make_media_record(items=[("small.png", 10), ("big.png", 200)])
        result = builder._ensure_media_size(record)
        assert result is not None
        assert len(result.value.items) == 1
        assert result.value.items[0].filename == "small.png"

    def test_all_items_exceed_size_returns_none(self, builder):
        """所有项都超过大小限制，返回 None"""
        record = _make_media_record(items=[("big1.png", 200), ("big2.png", 300)])
        result = builder._ensure_media_size(record)
        assert result is None

    def test_items_truncated_to_max_length(self, builder):
        """超过长度限制的项被截断"""
        record = _make_media_record(items=[("a.png", 10), ("b.png", 10), ("c.png", 10), ("d.png", 10)])
        result = builder._ensure_media_size(record)
        assert result is not None
        assert len(result.value.items) == 3
        assert [item.filename for item in result.value.items] == ["a.png", "b.png", "c.png"]

    def test_size_and_length_limits_combined(self, builder):
        """大小和长度限制同时生效：先过滤超大小，再截断超长度"""
        record = _make_media_record(
            items=[("a.png", 10), ("big.png", 200), ("b.png", 10), ("c.png", 10), ("d.png", 10)]
        )
        result = builder._ensure_media_size(record)
        # big.png 被过滤，剩余 a/b/c/d 共 4 项，截断为 3 项
        assert result is not None
        assert len(result.value.items) == 3
        assert [item.filename for item in result.value.items] == ["a.png", "b.png", "c.png"]

    def test_empty_record_returns_none(self, builder):
        """空 record（无 items）返回 None"""
        record = _make_media_record()
        result = builder._ensure_media_size(record)
        assert result is None

    def test_preserved_fields_in_new_record(self, builder):
        """新 record 保留了 key、step、type、timestamp 等字段"""
        record = _make_media_record(key="train/img", step=42, items=[("a.png", 10), ("big.png", 200)])
        from google.protobuf.timestamp_pb2 import Timestamp

        ts = Timestamp(seconds=1000)
        record.timestamp.CopyFrom(ts)
        result = builder._ensure_media_size(record)
        assert result is not None
        assert result.key == "train/img"
        assert result.step == 42
        assert result.timestamp == ts

    def test_exact_boundary_size_not_dropped(self, builder):
        """大小恰好等于限制的项不被丢弃"""
        record = _make_media_record(items=[("exact.png", 100)])
        result = builder._ensure_media_size(record)
        assert result is record

    def test_exact_boundary_length_not_truncated(self, builder):
        """大小恰好等于限制的项不被截断"""
        record = _make_media_record(items=[("a.png", 10), ("b.png", 10), ("c.png", 10)])
        result = builder._ensure_media_size(record)
        assert result is record


class TestResolveMediaDir:
    """_resolve_media_dir 依据 core.skip_store 决定媒体是否落盘"""

    def test_skip_store_returns_none_and_never_mkdir(self, tmp_path):
        """skip_store 下返回 None 且不创建 media 目录"""
        media_dir = tmp_path / "media"
        settings = SimpleNamespace(core=SimpleNamespace(skip_store=True))
        ctx = SimpleNamespace(config=SimpleNamespace(settings=settings), media_dir=media_dir)
        builder = RecordBuilder(ctx)  # type: ignore[arg-type]
        assert builder._resolve_media_dir(Text.column_type()) is None
        assert not media_dir.exists()


class TestBuildConfig:
    """build_config 依据 core.skip_store 决定内容内联进 payload 还是回读磁盘"""

    CONFIG_PATH = Path("/tmp/run/files/config.yaml")

    @staticmethod
    def _builder(skip_store: bool) -> RecordBuilder:
        settings = SimpleNamespace(core=SimpleNamespace(skip_store=skip_store))
        ctx = SimpleNamespace(config=SimpleNamespace(settings=settings))
        return RecordBuilder(ctx)  # type: ignore[arg-type]

    def _event(self) -> ConfigEvent:
        ts = Timestamp()
        ts.GetCurrentTime()
        return ConfigEvent(
            path=self.CONFIG_PATH,
            timestamp=ts,
            content={"lr": {"value": 0.01, "desc": "", "sort": 0}},
        )

    def test_skip_store_inlines_payload(self):
        """skip_store 下内容按落盘同款 YAML 编码填入 payload，source_path 留空"""
        event = self._event()
        record = self._builder(True).build_config(event)

        assert record.name == "config"
        assert record.type == SaveType.SAVE_TYPE_CONFIG
        assert record.source_path == ""
        assert yaml.safe_load(record.payload) == event.content

    def test_default_reads_from_disk(self):
        """默认模式 payload 恒空，由 Core 按 source_path 回读 config.yaml"""
        record = self._builder(False).build_config(self._event())

        assert record.source_path == self.CONFIG_PATH.as_posix()
        assert record.payload == b""


class TestIsScalarValue:
    """is_scalar_value 必须与 build_scalar_or_media 的分派结果一致。

    消费端用它做预扫描：若新增的媒体注册类型未同步 _NON_SCALAR_TYPES，
    预扫描会把媒体值当标量 transform，导致 media 预处理异常或行为漂移。
    """

    def test_registry_types_all_covered(self):
        """注册表中每个显式类型都必须在 _NON_SCALAR_TYPES 中"""
        # Python 3.12 下经 __get__ 包装后的方法不含 registry，需取类上的原始 descriptor
        sdm = RecordBuilder.__dict__["build_scalar_or_media"]
        registered = {tp for tp in sdm.dispatcher.registry if tp is not object}
        assert registered <= set(_NON_SCALAR_TYPES)

    def test_scalars_return_true(self):
        assert is_scalar_value(1.5) is True
        assert is_scalar_value("str") is True

    def test_media_returns_false(self):
        assert is_scalar_value([Text("a"), Text("b")]) is False
        assert is_scalar_value(Text("a")) is False

    def test_echarts_returns_false(self):
        """默认分支中 echarts 类型会被包装为 ECharts 媒体，非标量"""
        from pyecharts.charts import Line

        assert is_scalar_value(Line()) is False
