"""RetryPolicy 单元测试：错误分类、指数退避与 Retry-After 解析。"""

import email.utils
import math
from typing import Dict, Optional
from unittest.mock import MagicMock

import pytest

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy
from swanlab.sdk.typings.core_python.transport import FailureSource, UploadFailure


class _FakeApiErrorResponse:
    """ApiError 测试替身。"""

    def __init__(self, status_code: int, headers: Optional[Dict[str, str]] = None) -> None:
        self.status_code = status_code
        self.headers = headers or {}
        self.request = MagicMock()
        self.url = "https://api.test/upload"


def _api_failure(
    status: int, source: FailureSource = "swanlab_api", headers: Optional[Dict[str, str]] = None
) -> UploadFailure:
    error = ApiError(
        _FakeApiErrorResponse(status, headers),
        method="POST",
        trace_id="trace-id",
        code="test-error",
        message=f"http {status}",
    )
    return UploadFailure(source=source, message=f"http {status}", status=status, exception=error)


def _plain_failure(source: FailureSource, status: Optional[int] = None) -> UploadFailure:
    return UploadFailure(source=source, message="ConnectionError", status=status)


def test_classify_transient():
    """网络异常、5xx 与对象存储非 4xx 归类为 transient。"""
    policy = RetryPolicy()
    assert policy.classify(_plain_failure("swanlab_api")) == "transient"
    assert policy.classify(_api_failure(500)) == "transient"
    assert policy.classify(_api_failure(503)) == "transient"
    assert policy.classify(_api_failure(502, source="object_storage")) == "transient"


def test_classify_blocked():
    """4xx 拒绝、本地数据错误与 worker 内部错误归类为 blocked。"""
    policy = RetryPolicy()
    assert policy.classify(_api_failure(401)) == "blocked"
    assert policy.classify(_api_failure(403)) == "blocked"
    assert policy.classify(_api_failure(400)) == "blocked"
    assert policy.classify(_api_failure(404)) == "blocked"
    assert policy.classify(_plain_failure("local_data")) == "blocked"
    assert policy.classify(_plain_failure("worker_internal")) == "blocked"
    assert policy.classify(_api_failure(403, source="object_storage")) == "blocked"


def test_classify_rate_limited():
    assert RetryPolicy().classify(_api_failure(429)) == "rate_limited"


def test_backoff_grows_exponentially_and_caps():
    """指数退避按 2 倍增长并受上限 30 秒约束。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    delays = [policy.schedule([_plain_failure("swanlab_api")], round_).delay for round_ in range(1, 8)]
    assert delays == [1.0, 2.0, 4.0, 8.0, 16.0, 30.0, 30.0]


@pytest.mark.parametrize("failed_rounds", [1025, 1_000_000])
def test_backoff_does_not_overflow_on_huge_round_counts(failed_rounds):
    """大失败轮数下退避计算不溢出，结果不超过 max_delay。"""
    policy = RetryPolicy(rng=lambda _low, high_: high_)
    delay = policy.schedule([_plain_failure("swanlab_api")], failed_rounds).delay
    assert math.isfinite(delay)
    assert 0 < delay <= policy.max_delay


def test_backoff_jitter_within_bounds_and_restricted_by_cap():
    """随机抖动限制在 0.8~1.2 倍且不超过 max_delay。"""
    low = RetryPolicy(rng=lambda low_, _high: low_)
    high = RetryPolicy(rng=lambda _low, high_: high_)
    assert low.schedule([_plain_failure("swanlab_api")], 1).delay == 0.8
    assert high.schedule([_plain_failure("swanlab_api")], 1).delay == 1.2
    assert high.schedule([_plain_failure("swanlab_api")], 6).delay == 30.0


def test_blocked_delay_is_flat_and_ignores_round_count():
    """blocked 间隔固定为 60 秒。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.2)
    for round_ in (1, 3, 10):
        assert policy.schedule([_api_failure(401)], round_).delay == 60.0
        assert policy.schedule([_api_failure(401)], round_).kind == "blocked"


def test_schedule_empty_failures_is_blocked():
    """空失败列表按 blocked 处理。"""
    schedule = RetryPolicy().schedule([], 1)
    assert schedule.kind == "blocked"
    assert schedule.delay == 60.0


def test_rate_limited_uses_retry_after_seconds_uncapped():
    """429 秒数 Retry-After 不受普通上限截短，并回填 failure.retry_after。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    failure = _api_failure(429, headers={"Retry-After": "120"})
    schedule = policy.schedule([failure], 1)
    assert schedule.kind == "rate_limited"
    assert schedule.delay == 120.0
    assert failure.retry_after == 120.0


def test_rate_limited_http_date_converted_via_wall_clock():
    """429 HTTP 日期格式通过 wall_clock 换算为等待秒数。"""
    wall_now = 1_700_000_000.0
    policy = RetryPolicy(wall_clock=lambda: wall_now, rng=lambda _low, _high: 1.0)
    header = email.utils.formatdate(wall_now + 90, usegmt=True)
    schedule = policy.schedule([_api_failure(429, headers={"Retry-After": header})], 1)
    assert schedule.delay == 90.0


def test_rate_limited_invalid_or_missing_header_falls_back_to_backoff():
    """缺失或无效的 Retry-After 回退至指数退避。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    for headers in (None, {"Retry-After": "not-a-date"}, {"Retry-After": "-5"}):
        failure = _api_failure(429, headers=headers)
        schedule = policy.schedule([failure], 2)
        assert schedule.kind == "transient"
        assert schedule.delay == 2.0


@pytest.mark.parametrize("value", [None, "", "   ", "garbage"])
def test_parse_retry_after_invalid_forms(value):
    assert RetryPolicy().parse_retry_after(value) is None


def test_parse_retry_after_zero_is_valid():
    assert RetryPolicy().parse_retry_after("0") == 0.0


def test_http_date_in_past_is_invalid():
    """过去的 HTTP 日期视为无效。"""
    wall_now = 1_700_000_000.0
    policy = RetryPolicy(wall_clock=lambda: wall_now)
    header = email.utils.formatdate(wall_now - 60, usegmt=True)
    assert policy.parse_retry_after(header) is None


def test_parse_retry_after_rejects_non_finite():
    """非有限数值（nan/inf/1e999）视为无效。"""
    policy = RetryPolicy()
    for value in ("nan", "inf", "-inf", "Infinity", "1e999"):
        assert policy.parse_retry_after(value) is None


def test_rate_limited_non_finite_header_falls_back_to_backoff():
    """429 非有限 Retry-After 回退至指数退避。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    for value in ("nan", "inf", "1e999"):
        failure = _api_failure(429, headers={"Retry-After": value})
        schedule = policy.schedule([failure], 2)
        assert schedule.kind == "transient"
        assert schedule.delay == 2.0


def test_rate_limited_huge_finite_retry_after_is_honored():
    """大有限值 Retry-After 原样保留。"""
    failure = _api_failure(429, headers={"Retry-After": "1e20"})
    schedule = RetryPolicy().schedule([failure], 1)
    assert schedule.kind == "rate_limited"
    assert schedule.delay == 1e20


def test_multiple_failures_take_longest_constraint():
    """多个失败取最大等待时长。"""
    policy = RetryPolicy(rng=lambda low_, _high: low_)
    transient = _api_failure(503)
    limited = _api_failure(429, headers={"Retry-After": "45"})
    schedule = policy.schedule([transient, limited], 1)
    assert schedule.delay == 45.0


def test_retry_after_takes_precedence_over_backoff():
    """有效 Retry-After 不与指数退避混合计算。"""
    policy = RetryPolicy(rng=lambda _low, high_: high_)
    limited = _api_failure(429, headers={"Retry-After": "3"})
    assert policy.schedule([limited], 1).delay == 3.0
    assert policy.schedule([limited], 5).delay == 3.0
