"""
@author: nexisato
@file: test_sync_range.py
@time: 2026/10/23
@description: sync 范围过滤器与时间解析测试
"""

from datetime import datetime, timedelta, timezone

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import SyncRange
from swanlab.sdk.internal.core_python.sync_range import SyncRangeFilter
from swanlab.utils.time import parse_timestamp_ns, parse_timezone


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


def test_parse_timezone_supports_utc_and_offsets():
    assert parse_timezone(None) is None
    assert parse_timezone("UTC") is timezone.utc
    assert parse_timezone("Z") is timezone.utc
    assert parse_timezone("+08:00") == timezone(timedelta(hours=8))
    assert parse_timezone("-0530") == timezone(timedelta(hours=-5, minutes=-30))
    assert parse_timezone("+8") == timezone(timedelta(hours=8))


def test_parse_timezone_rejects_invalid_values():
    for value in ("", "+25:00", "+08:70", "not/a-timezone"):
        with pytest.raises(ValueError):
            parse_timezone(value)


def test_parse_timezone_supports_iana_names():
    try:
        parsed = parse_timezone("Asia/Shanghai")
    except ValueError:
        pytest.skip("system timezone database is unavailable")
    assert parsed is not None
    assert parsed.utcoffset(datetime(2024, 8, 1)) == timedelta(hours=8)


def test_parse_timestamp_ns_applies_given_timezone():
    expected = 1722441600 * 1_000_000_000  # 2024-08-01T00:00:00+08:00

    assert parse_timestamp_ns("2024-08-01 00:00:00", tz=timezone(timedelta(hours=8))) == expected
    assert parse_timestamp_ns(datetime(2024, 8, 1, 0, 0, 0), tz=timezone(timedelta(hours=8))) == expected


def test_parse_timestamp_ns_keeps_explicit_offset():
    expected = 1722470400 * 1_000_000_000  # 2024-08-01T00:00:00Z

    # 输入自带时区时，tz 参数不覆盖输入
    assert parse_timestamp_ns("2024-08-01T00:00:00Z", tz=timezone(timedelta(hours=8))) == expected
