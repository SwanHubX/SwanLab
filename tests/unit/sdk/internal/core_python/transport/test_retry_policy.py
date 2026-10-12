"""RetryPolicy 单元测试：指数退避。"""

import math

import pytest

from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy
from swanlab.sdk.typings.core_python.transport import FailureSource, UploadFailure


def _failure(source: FailureSource = "swanlab_api") -> UploadFailure:
    return UploadFailure(source=source, message="error")


def test_backoff_grows_exponentially_and_caps():
    """指数退避按 2 倍增长并受上限 30 秒约束。"""
    policy = RetryPolicy(rng=lambda low_, high_: 1.0)
    delays = [policy.schedule([_failure()], round_).delay for round_ in range(1, 8)]
    assert delays == [1.0, 2.0, 4.0, 8.0, 16.0, 30.0, 30.0]


@pytest.mark.parametrize("failed_rounds", [1025, 1_000_000])
def test_backoff_does_not_overflow_on_huge_round_counts(failed_rounds):
    """大失败轮数下退避计算不溢出，结果不超过 max_delay。"""
    policy = RetryPolicy(rng=lambda low_, high_: high_)
    delay = policy.schedule([_failure()], failed_rounds).delay
    assert math.isfinite(delay)
    assert 0 < delay <= policy.max_delay


def test_backoff_jitter_within_bounds_and_restricted_by_cap():
    """随机抖动限制在 0.8~1.2 倍且不超过 max_delay。"""
    low = RetryPolicy(rng=lambda low_, high_: low_)
    high = RetryPolicy(rng=lambda low_, high_: high_)
    assert low.schedule([_failure()], 1).delay == 0.8
    assert high.schedule([_failure()], 1).delay == 1.2
    assert high.schedule([_failure()], 6).delay == 30.0


def test_schedule_empty_failures_uses_base_delay():
    """空失败列表按第 1 轮退避处理。"""
    policy = RetryPolicy(rng=lambda low_, high_: 1.0)
    assert policy.schedule([], 1).delay == 1.0


def test_schedule_ignores_failure_source():
    """所有失败来源（包括 local_data）均走相同的指数退避。"""
    policy = RetryPolicy(rng=lambda low_, high_: 1.0)
    for source in ("swanlab_api", "object_storage", "local_data", "worker_internal"):
        assert policy.schedule([_failure(source)], 1).delay == 1.0  # type: ignore[arg-type]
