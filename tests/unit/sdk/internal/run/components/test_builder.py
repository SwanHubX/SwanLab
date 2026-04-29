"""
@author: cunyue
@file: test_builder.py
@time: 2026/4/29
@description: RecordBuilder 单元测试
"""

import pytest

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem, MediaRecord
from swanlab.sdk.internal.run.components.builder import RecordBuilder


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
        """长度恰好等于限制的项不被截断"""
        record = _make_media_record(items=[("a.png", 10), ("b.png", 10), ("c.png", 10)])
        result = builder._ensure_media_size(record)
        assert result is record
