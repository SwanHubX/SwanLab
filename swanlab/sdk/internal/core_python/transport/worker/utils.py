"""
@author: caddiesnew
@file: utils.py
@time: 2026/10/11
@description: worker 重试策略：错误分类、指数退避、Retry-After 解析
"""

import math
import random
import time
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional, Sequence, Tuple

from swanlab.exceptions import ApiError
from swanlab.sdk.typings.core_python.transport import UploadFailure

RetryKind = Literal["transient", "rate_limited", "blocked"]


@dataclass(frozen=True)
class RetrySchedule:
    """单轮失败后的重试调度结果。"""

    kind: RetryKind
    delay: float


class RetryPolicy:
    """根据上传失败信息计算下次重试等待时长。"""

    BASE_DELAY: float = 1.0
    MULTIPLIER: float = 2.0
    MAX_DELAY: float = 30.0
    JITTER_RANGE: Tuple[float, float] = (0.8, 1.2)
    BLOCKED_DELAY: float = 60.0

    def __init__(
        self,
        *,
        base_delay: Optional[float] = None,
        multiplier: Optional[float] = None,
        max_delay: Optional[float] = None,
        blocked_delay: Optional[float] = None,
        clock: Callable[[], float] = time.monotonic,
        wall_clock: Callable[[], float] = time.time,
        rng: Optional[Callable[[float, float], float]] = None,
    ) -> None:
        self.base_delay = base_delay if base_delay is not None else self.BASE_DELAY
        self.multiplier = multiplier if multiplier is not None else self.MULTIPLIER
        self.max_delay = max_delay if max_delay is not None else self.MAX_DELAY
        self.blocked_delay = blocked_delay if blocked_delay is not None else self.BLOCKED_DELAY
        self.clock = clock
        self.wall_clock = wall_clock
        self._rng = rng if rng is not None else random.uniform
        # 指数上限，防止失败轮数过大时幂运算溢出
        if self.multiplier > 1.0 and self.max_delay > self.base_delay:
            self._saturation_exponent = math.ceil(math.log(self.max_delay / self.base_delay, self.multiplier))
        else:
            self._saturation_exponent = 0

    def classify(self, failure: UploadFailure) -> RetryKind:
        """对单个上传失败进行重试分类。"""
        if failure.source in ("local_data", "worker_internal"):
            return "blocked"
        if failure.source == "object_storage":
            if failure.status is not None and 400 <= failure.status < 500:
                return "blocked"
            return "transient"
        if failure.status == 429:
            return "rate_limited"
        if failure.status is not None and 400 <= failure.status < 500:
            return "blocked"
        return "transient"

    def schedule(self, failures: Sequence[UploadFailure], failed_rounds: int) -> RetrySchedule:
        """根据本轮失败列表与连续失败轮数计算重试间隔，多个失败取最大等待时长。"""
        if not failures:
            return RetrySchedule("blocked", self.blocked_delay)
        schedules: List[RetrySchedule] = []
        for failure in failures:
            kind = self.classify(failure)
            if kind == "rate_limited":
                retry_after = self.retry_after_of(failure)
                if retry_after is not None:
                    failure.retry_after = retry_after
                    schedules.append(RetrySchedule(kind, retry_after))
                    continue
                kind = "transient"
            delay = self._backoff_delay(failed_rounds) if kind == "transient" else self.blocked_delay
            schedules.append(RetrySchedule(kind, delay))
        return max(schedules, key=lambda item: item.delay)

    def _backoff_delay(self, failed_rounds: int) -> float:
        """计算带随机抖动的指数退避等待时长。"""
        exponent = max(failed_rounds, 1) - 1
        exponent = min(exponent, self._saturation_exponent)
        base = min(self.base_delay * (self.multiplier**exponent), self.max_delay)
        low, high = self.JITTER_RANGE
        return min(base * self._rng(low, high), self.max_delay)

    def retry_after_of(self, failure: UploadFailure) -> Optional[float]:
        """从 ApiError 响应头解析 Retry-After。"""
        if not isinstance(failure.exception, ApiError):
            return None
        response = getattr(failure.exception, "response", None)
        header = getattr(response, "headers", {}).get("Retry-After")
        return self.parse_retry_after(header)

    def parse_retry_after(self, value: Optional[str]) -> Optional[float]:
        """将 Retry-After 字符串解析为等待秒数（支持秒数与 HTTP 日期格式）。"""
        if value is None or not value.strip():
            return None
        text = value.strip()
        try:
            seconds = float(text)
        except ValueError:
            date = _parse_http_date(text)
            if date is None:
                return None
            seconds = date - self.wall_clock()
        if seconds < 0 or not math.isfinite(seconds):
            return None
        return seconds


def _parse_http_date(text: str) -> Optional[float]:
    """将 RFC 7231 HTTP 日期字符串解析为 Unix 时间戳。"""
    import email.utils
    from datetime import timezone

    try:
        date = email.utils.parsedate_to_datetime(text)
    except (TypeError, ValueError):
        return None
    if date is None:
        return None
    if date.tzinfo is None:
        date = date.replace(tzinfo=timezone.utc)
    return date.timestamp()


__all__ = [
    "RetryKind",
    "RetryPolicy",
    "RetrySchedule",
]
