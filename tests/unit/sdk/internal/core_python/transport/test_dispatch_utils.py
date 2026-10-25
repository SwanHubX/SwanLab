"""Dispatch 分片与分组工具测试。"""

import pytest

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.dispatch import generate_chunks, group_records_by_type


def test_generate_chunks_single_chunk(make_scalar_record):
    """记录数不超过 batch_size 时产出单个分片。"""
    records = [make_scalar_record(step=i) for i in range(5)]
    chunks = list(generate_chunks(records, batch_size=10))

    assert len(chunks) == 1
    assert chunks[0] == (records, 5)


def test_generate_chunks_exact_boundary(make_scalar_record):
    """记录数等于 batch_size 时产出单个分片。"""
    records = [make_scalar_record(step=i) for i in range(100)]
    chunks = list(generate_chunks(records, batch_size=100))

    assert len(chunks) == 1
    assert chunks[0] == (records, 100)


def test_generate_chunks_multiple_chunks(make_scalar_record):
    """记录数超过 batch_size 时拆分为多个分片。"""
    records = [make_scalar_record(step=i) for i in range(25)]
    chunks = list(generate_chunks(records, batch_size=10))

    assert len(chunks) == 3
    assert chunks[0][1] == 10
    assert chunks[1][1] == 10
    assert chunks[2][1] == 5


def test_generate_chunks_no_limit(make_scalar_record):
    """batch_size=-1 时返回单个分片。"""
    records = [make_scalar_record(step=i) for i in range(9999)]
    chunks = list(generate_chunks(records, batch_size=-1))

    assert len(chunks) == 1
    assert chunks[0] == (records, 9999)


def test_generate_chunks_invalid_max_records_per_request(make_scalar_record):
    """batch_size 为 0 或小于 -1 时抛出 ValueError。"""
    records = [make_scalar_record(step=1)]

    with pytest.raises(ValueError, match="max_records_per_request"):
        list(generate_chunks(records, batch_size=0))

    with pytest.raises(ValueError, match="max_records_per_request"):
        list(generate_chunks(records, batch_size=-5))


def test_generate_chunks_empty_records():
    """空记录列表产出空迭代器。"""
    chunks = list(generate_chunks([], batch_size=10))
    assert len(chunks) == 0


def test_generate_chunks_no_limit_empty():
    """batch_size=-1 且列表为空时产出空迭代器。"""
    chunks = list(generate_chunks([], batch_size=-1))
    assert len(chunks) == 0


def test_group_records_by_type_basic(make_scalar_record, make_config_record):
    """按 record_type 对 Record 分组。"""
    metric_1 = make_scalar_record(step=1)
    config = make_config_record()
    metric_2 = make_scalar_record(step=2)

    grouped, invalid = group_records_by_type([metric_1, config, metric_2])

    assert invalid == []
    assert set(grouped.keys()) == {"scalar", "save"}
    assert grouped["scalar"] == [metric_1, metric_2]
    assert grouped["save"] == [config]


def test_group_records_by_type_preserves_order(make_scalar_record, make_config_record):
    """分组 key 保持首次出现顺序。"""
    config = make_config_record()
    metric = make_scalar_record(step=1)

    grouped, _ = group_records_by_type([config, metric])
    assert list(grouped.keys()) == ["save", "scalar"]


def test_group_records_by_type_empty():
    """空列表返回空字典。"""
    assert group_records_by_type([]) == ({}, [])


def test_group_records_by_type_rejects_non_record():
    """传入非 Record 实例时抛出 TypeError。"""
    with pytest.raises(TypeError, match="Record"):
        group_records_by_type(["not_a_record"])  # type: ignore[list-item]


def test_group_records_by_type_collects_untypeable_record():
    """未设置 record_type 的空 Record 被单独返回，不影响有效分组。"""
    empty = Record()
    grouped, invalid = group_records_by_type([empty])

    assert grouped == {}
    assert invalid == [empty]
