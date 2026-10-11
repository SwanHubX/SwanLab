"""RetryPolicy 契约测试：错误分类、指数退避与抖动、Retry-After 解析（plan §6）。

全部通过注入 rng / wall_clock 完成，不依赖真实 sleep；worker 侧的等待与唤醒行为
见 test_transport.py。
"""

import email.utils
from typing import Dict, Optional
from unittest.mock import MagicMock

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy
from swanlab.sdk.typings.core_python.transport import FailureSource, UploadFailure


class _FakeApiErrorResponse:
    """构造 ApiError 所需的最小 Response 替身，带可配置响应头。"""

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
    """无 ApiError 异常的失败（网络异常等）。"""
    return UploadFailure(source=source, message="ConnectionError", status=status)


# ── 分类（plan §6.1） ────────────────────────────────────────


def test_classify_transient():
    """网络异常（无状态码）与 5xx、对象存储非 4xx 属可重试。"""
    policy = RetryPolicy()
    assert policy.classify(_plain_failure("swanlab_api")) == "transient"
    assert policy.classify(_api_failure(500)) == "transient"
    assert policy.classify(_api_failure(503)) == "transient"
    assert policy.classify(_api_failure(502, source="object_storage")) == "transient"


def test_classify_blocked():
    """401/403 与其他业务拒绝、本地坏数据、worker 内部错误、对象存储 4xx 属低频受阻。"""
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


# ── 指数退避（plan §6.2） ────────────────────────────────────


def test_backoff_grows_exponentially_and_caps():
    """抖动固定为 1.0 时：1、2、4、8、16、30、30——2 倍增长，上限 30 秒。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    delays = [policy.schedule([_plain_failure("swanlab_api")], round_).delay for round_ in range(1, 8)]
    assert delays == [1.0, 2.0, 4.0, 8.0, 16.0, 30.0, 30.0]


def test_backoff_jitter_within_bounds_and_restricted_by_cap():
    """抖动为计算值的 0.8～1.2 倍，且乘抖动后再限制到普通上限。"""
    low = RetryPolicy(rng=lambda low_, _high: low_)
    high = RetryPolicy(rng=lambda _low, high_: high_)
    # 第 1 轮 base=1：抖动下界 0.8，上界 1.2
    assert low.schedule([_plain_failure("swanlab_api")], 1).delay == 0.8
    assert high.schedule([_plain_failure("swanlab_api")], 1).delay == 1.2
    # 第 6 轮 base=30（已到上限）：乘 1.2 后仍被限制回 30
    assert high.schedule([_plain_failure("swanlab_api")], 6).delay == 30.0


def test_blocked_delay_is_flat_and_ignores_round_count():
    """低频受阻间隔固定 60 秒，不随连续失败轮数增长、不受抖动影响。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.2)
    for round_ in (1, 3, 10):
        assert policy.schedule([_api_failure(401)], round_).delay == 60.0
        assert policy.schedule([_api_failure(401)], round_).kind == "blocked"


def test_schedule_empty_failures_is_blocked():
    """worker 兜底路径（失败原因为空）：按低频受阻处理，避免空转。"""
    schedule = RetryPolicy().schedule([], 1)
    assert schedule.kind == "blocked"
    assert schedule.delay == 60.0


# ── Retry-After（plan §6.2：有效值不被普通上限截短） ─────────


def test_rate_limited_uses_retry_after_seconds_uncapped():
    """429 携带秒数 Retry-After：直接采用，120 秒不被 30 秒上限截短，且回填 failure。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    failure = _api_failure(429, headers={"Retry-After": "120"})
    schedule = policy.schedule([failure], 1)
    assert schedule.kind == "rate_limited"
    assert schedule.delay == 120.0
    assert failure.retry_after == 120.0


def test_rate_limited_http_date_converted_via_wall_clock():
    """429 携带 HTTP 日期：墙钟换算一次为相对时长。"""
    wall_now = 1_700_000_000.0
    policy = RetryPolicy(wall_clock=lambda: wall_now, rng=lambda _low, _high: 1.0)
    header = email.utils.formatdate(wall_now + 90, usegmt=True)
    schedule = policy.schedule([_api_failure(429, headers={"Retry-After": header})], 1)
    assert schedule.delay == 90.0


def test_rate_limited_invalid_or_missing_header_falls_back_to_backoff():
    """缺失、畸形、负值 Retry-After 回退普通退避。"""
    policy = RetryPolicy(rng=lambda _low, _high: 1.0)
    for headers in (None, {"Retry-After": "not-a-date"}, {"Retry-After": "-5"}):
        failure = _api_failure(429, headers=headers)
        schedule = policy.schedule([failure], 2)
        assert schedule.kind == "transient"
        assert schedule.delay == 2.0


def test_parse_retry_after_invalid_forms():
    policy = RetryPolicy()
    assert policy.parse_retry_after(None) is None
    assert policy.parse_retry_after("") is None
    assert policy.parse_retry_after("   ") is None
    assert policy.parse_retry_after("garbage") is None


def test_parse_retry_after_zero_is_valid():
    """Retry-After: 0 合法（服务端允许立即重试）。"""
    assert RetryPolicy().parse_retry_after("0") == 0.0


def test_http_date_in_past_is_invalid():
    """过去的 HTTP 日期换算为负值，视为无效，回退普通退避。"""
    wall_now = 1_700_000_000.0
    policy = RetryPolicy(wall_clock=lambda: wall_now)
    header = email.utils.formatdate(wall_now - 60, usegmt=True)
    assert policy.parse_retry_after(header) is None


# ── 多失败取最长约束（plan §6.2） ────────────────────────────


def test_multiple_failures_take_longest_constraint():
    """多个失败取最长约束；随机抖动不得使结果低于有效 Retry-After。"""
    policy = RetryPolicy(rng=lambda low_, _high: low_)  # 抖动取下界 0.8
    transient = _api_failure(503)
    limited = _api_failure(429, headers={"Retry-After": "45"})
    schedule = policy.schedule([transient, limited], 1)
    assert schedule.delay == 45.0


def test_retry_after_takes_precedence_over_backoff():
    """有效 Retry-After 优先于指数退避（plan §6.1）：不与退避取 max，也不被抖动改写。"""
    policy = RetryPolicy(rng=lambda _low, high_: high_)  # 抖动取上界 1.2
    limited = _api_failure(429, headers={"Retry-After": "3"})
    # 无论退避轮数多大（第 5 轮退避 16*1.2=19.2），有效 Retry-After 独占生效
    assert policy.schedule([limited], 1).delay == 3.0
    assert policy.schedule([limited], 5).delay == 3.0
