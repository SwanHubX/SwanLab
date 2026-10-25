"""RetryPolicy 单元测试：指数退避。"""

import math

from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy


def test_backoff_grows_exponentially_and_caps():
    """指数退避按 2 倍增长并受上限 10 秒约束。"""
    policy = RetryPolicy(rng=lambda low_, high_: 1.0)
    delays = [policy.schedule(round_) for round_ in range(1, 8)]
    assert delays == [1.0, 2.0, 4.0, 8.0, 10.0, 10.0, 10.0]


def test_backoff_does_not_overflow_on_huge_round_counts():
    """大失败轮数下退避计算不溢出，结果不超过 max_delay。"""
    policy = RetryPolicy(rng=lambda low_, high_: high_)
    delay = policy.schedule(1_000_000)
    assert math.isfinite(delay)
    assert 0 < delay <= policy.max_delay


def test_backoff_jitter_within_bounds_and_restricted_by_cap():
    """随机抖动限制在 0.8~1.2 倍且不超过 max_delay。"""
    low = RetryPolicy(rng=lambda low_, high_: low_)
    high = RetryPolicy(rng=lambda low_, high_: high_)
    assert low.schedule(1) == 0.8
    assert high.schedule(1) == 1.2
    assert high.schedule(6) == 10.0
