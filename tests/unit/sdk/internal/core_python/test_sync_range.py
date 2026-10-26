"""
@author: nexisato
@file: test_sync_range.py
@time: 2026/10/23
@description: sync 范围过滤器与时间解析测试
"""

from datetime import datetime, timezone

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import SyncRange
from swanlab.sdk.internal.core_python.sync_range import SyncRangeFilter
from swanlab.utils.time import parse_timestamp_ns


def make_timestamp(seconds: int, nanos: int = 0) -> Timestamp:
    ts = Timestamp()
    ts.seconds = seconds
    ts.nanos = nanos
    return ts


def test_filter_disabled_without_range():
    f = SyncRangeFilter.from_proto(None)

    assert f.enabled is False
    assert f.match(step=1, timestamp=make_timestamp(1)) is True


def test_filter_disabled_when_timestamps_are_zero():
    f = SyncRangeFilter.from_proto(SyncRange(start_time=Timestamp(), end_time=Timestamp()))

    assert f.enabled is False


def test_filter_step_range_is_inclusive():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=10, end_step=20))

    assert f.enabled is True
    assert f.match(step=10) is True
    assert f.match(step=20) is True
    assert f.match(step=9) is False
    assert f.match(step=21) is False


def test_filter_step_range_ignores_records_without_step():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=10, end_step=20))

    # 终端日志没有 step，只受时间范围约束
    assert f.match(timestamp=make_timestamp(1)) is True


def test_filter_reports_time_bounds():
    assert SyncRangeFilter.from_proto(SyncRange(start_step=1)).has_time is False
    assert SyncRangeFilter.from_proto(SyncRange(start_time=make_timestamp(1))).has_time is True


def test_filter_open_ended_step_range():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=100))

    assert f.match(step=100) is True
    assert f.match(step=10_000) is True
    assert f.match(step=99) is False


def test_filter_time_range_is_inclusive():
    f = SyncRangeFilter.from_proto(SyncRange(start_time=make_timestamp(100, 500), end_time=make_timestamp(200)))

    assert f.match(step=1, timestamp=make_timestamp(100, 500)) is True
    assert f.match(step=1, timestamp=make_timestamp(200)) is True
    assert f.match(step=1, timestamp=make_timestamp(100, 499)) is False
    assert f.match(step=1, timestamp=make_timestamp(201)) is False


def test_filter_combines_step_and_time():
    f = SyncRangeFilter.from_proto(SyncRange(start_step=2, end_step=3, start_time=make_timestamp(100)))

    assert f.match(step=2, timestamp=make_timestamp(100)) is True
    assert f.match(step=1, timestamp=make_timestamp(100)) is False
    assert f.match(step=2, timestamp=make_timestamp(99)) is False


def test_parse_timestamp_ns_accepts_iso_and_unix():
    expected = 1722470400 * 1_000_000_000

    assert parse_timestamp_ns("2024-08-01T00:00:00Z") == expected
    assert parse_timestamp_ns(1722470400) == expected
    assert parse_timestamp_ns(1722470400000) == expected
    assert parse_timestamp_ns("1722470400") == expected


def test_parse_timestamp_ns_treats_naive_datetime_as_local():
    naive = datetime(2024, 8, 1, 8, 30, 0)

    assert parse_timestamp_ns(naive) == int(round(naive.timestamp() * 1_000_000_000))
    assert parse_timestamp_ns("2024-08-01 08:30:00") == int(round(naive.timestamp() * 1_000_000_000))


@pytest.mark.parametrize("value", ["", "  ", "not-a-time", -1, 0, float("nan")])
def test_parse_timestamp_ns_rejects_invalid_values(value):
    with pytest.raises(ValueError):
        parse_timestamp_ns(value)


def test_parse_timestamp_ns_rejects_bool():
    with pytest.raises(TypeError):
        parse_timestamp_ns(True)


def test_parse_timestamp_ns_keeps_subsecond_precision():
    """亚秒输入必须精确到纳秒，否则闭区间端点会被判到范围外。"""
    base = 1722470400 * 1_000_000_000

    assert parse_timestamp_ns("2024-08-01T00:00:00.000001Z") == base + 1_000
    assert parse_timestamp_ns("2024-08-01T00:00:00.000002Z") == base + 2_000
    assert parse_timestamp_ns(datetime(2024, 8, 1, 0, 0, 0, 1, tzinfo=timezone.utc)) == base + 1_000
    assert parse_timestamp_ns(1722470400123) == base + 123_000_000
    assert parse_timestamp_ns("1722470400123") == base + 123_000_000
    assert parse_timestamp_ns(1722470400.001) == base + 1_000_000


def test_filter_keeps_subsecond_endpoint():
    """从同一字面量解析出的时间戳必须命中它自己构成的闭区间。"""
    instant = parse_timestamp_ns("2024-08-01T00:00:00.000001Z")
    timestamp = Timestamp()
    timestamp.FromNanoseconds(instant)

    f = SyncRangeFilter(start_ns=instant, end_ns=instant)

    assert f.match(timestamp=timestamp) is True
