import pytest

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.helper import generate_chunks, group_records_by_type

# ─────────────────── generate_chunks ───────────────────


def test_generate_chunks_single_chunk(make_scalar_record):
    """验证记录数不超过 per_request_len 时产出单个分片。"""
    records = [make_scalar_record(step=i) for i in range(5)]
    chunks = list(generate_chunks(records, per_request_len=10))

    assert len(chunks) == 1
    assert chunks[0] == (records, 5)

def test_generate_chunks_multiple_chunks(make_scalar_record):
    """验证记录数超过 per_request_len 时正确拆分为多个分片，最后一片长度不足也正常产出。"""
    records = [make_scalar_record(step=i) for i in range(25)]
    chunks = list(generate_chunks(records, per_request_len=10))

    assert len(chunks) == 3
    assert chunks[0][1] == 10
    assert chunks[1][1] == 10
    assert chunks[2][1] == 5


def test_generate_chunks_no_limit(make_scalar_record):
    """验证 per_request_len=-1 时不拆分，一次性返回全部记录。"""
    records = [make_scalar_record(step=i) for i in range(9999)]
    chunks = list(generate_chunks(records, per_request_len=-1))

    assert len(chunks) == 1
    assert chunks[0] == (records, 9999)


def test_generate_chunks_invalid_per_request_len(make_scalar_record):
    """验证 per_request_len 为 0 或负数（非 -1）时抛出 ValueError。"""
    records = [make_scalar_record(step=1)]

    with pytest.raises(ValueError, match="per_request_len"):
        list(generate_chunks(records, per_request_len=0))

    with pytest.raises(ValueError, match="per_request_len"):
        list(generate_chunks(records, per_request_len=-5))


def test_generate_chunks_empty_records():
    """验证空记录列表在 per_request_len > 0 时产出单个空分片。"""
    chunks = list(generate_chunks([], per_request_len=10))
    assert len(chunks) == 1
    assert chunks[0] == ([], 0)


def test_generate_chunks_no_limit_empty():
    """验证 per_request_len=-1 且空列表时产出单个空分片。"""
    chunks = list(generate_chunks([], per_request_len=-1))
    assert len(chunks) == 1
    assert chunks[0] == ([], 0)


# ─────────────────── group_records_by_type ───────────────────


def test_group_records_by_type_basic(make_scalar_record, make_config_record):
    """验证不同 record_type 的 Record 被正确分组到各自的列表中。"""
    metric_1 = make_scalar_record(step=1)
    config = make_config_record()
    metric_2 = make_scalar_record(step=2)

    grouped = group_records_by_type([metric_1, config, metric_2])

    assert set(grouped.keys()) == {"metric", "config"}
    assert grouped["metric"] == [metric_1, metric_2]
    assert grouped["config"] == [config]


def test_group_records_by_type_preserves_order(make_scalar_record, make_config_record):
    """验证分组的 key 按 Record 首次出现顺序排列（OrderedDict）。"""
    config = make_config_record()
    metric = make_scalar_record(step=1)

    grouped = group_records_by_type([config, metric])
    assert list(grouped.keys()) == ["config", "metric"]


def test_group_records_by_type_empty():
    """验证空列表返回空字典。"""
    assert group_records_by_type([]) == {}


def test_group_records_by_type_rejects_non_record():
    """验证传入非 Record 实例时抛出 TypeError。"""
    with pytest.raises(TypeError, match="Record"):
        group_records_by_type(["not_a_record"])  # type: ignore[list-item]


def test_group_records_by_type_rejects_empty_record():
    """验证未设置任何 record_type 的空 Record 抛出 ValueError。"""
    with pytest.raises(ValueError, match="no active record_type"):
        group_records_by_type([Record()])
