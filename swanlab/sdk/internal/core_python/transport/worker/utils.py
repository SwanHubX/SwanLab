"""
@author: caddiesnew
@file: utils.py
@time: 2026/10/11
@description: worker 重试策略：错误分类、指数退避、Retry-After 解析（纯计算，可注入时钟与随机源）

本模块只做"一轮失败 → 下次尝试间隔"的计算，不做等待、不持有状态；
等待与唤醒由 worker 的 Condition 循环负责（plan §6.3）。
"""

import random
import time
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional, Sequence, Tuple

from swanlab.exceptions import ApiError
from swanlab.sdk.typings.core_python.transport import UploadFailure

# 调度类别：transient 指数退避、rate_limited 用 Retry-After、blocked 固定低频间隔
RetryKind = Literal["transient", "rate_limited", "blocked"]


@dataclass(frozen=True)
class RetrySchedule:
    """一轮失败后的重试调度结果：kind 用于告警文案，delay 为下次尝试的等待秒数。"""

    kind: RetryKind
    delay: float


class RetryPolicy:
    """失败轮 → 下次尝试间隔的纯计算策略（plan §6）。

    - transient（网络异常、5xx、对象存储传输错误）：指数退避 + 抖动，普通上限 30s
    - rate_limited（429）：有效 Retry-After 优先且不被普通上限截短，否则按 transient 退避
    - blocked（认证/权限/业务拒绝、本地坏数据、worker 内部错误）：固定 60s 低频尝试，
      不推断为可丢弃；本地坏数据在文件恢复后仍可自愈

    时钟与随机源可注入（plan §6.2），测试不依赖真实长时间 sleep。
    """

    # 普通失败首次等待（秒）
    BASE_DELAY: float = 1.0
    # 指数倍率
    MULTIPLIER: float = 2.0
    # 普通失败等待上限（秒），不截短有效 Retry-After
    MAX_DELAY: float = 30.0
    # 抖动区间：计算值的该倍数范围内取随机，再限制到普通上限
    JITTER_RANGE: Tuple[float, float] = (0.8, 1.2)
    # 认证/业务拒绝/本地受阻的低频尝试间隔（秒）
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
        # 单调时钟：退避截止时间计算；墙钟：HTTP 日期换算（plan §6.2）
        self.clock = clock
        self.wall_clock = wall_clock
        self._rng = rng if rng is not None else random.uniform

    # ── 错误分类（plan §6.1） ──

    def classify(self, failure: UploadFailure) -> RetryKind:
        """单个失败 → 调度类别。403 不统一解释为 Key 撤销，保留来源与业务码由告警展示。"""
        if failure.source in ("local_data", "worker_internal"):
            # 本地坏数据/内部错误属确定性受阻：低频尝试避免空转，源文件恢复后可自愈
            return "blocked"
        if failure.source == "object_storage":
            # 对象存储 4xx 无法可靠区分签名失效与权限受限（地址在整流程重跑时天然重新签名），低频等待
            if failure.status is not None and 400 <= failure.status < 500:
                return "blocked"
            return "transient"
        if failure.status == 429:
            return "rate_limited"
        if failure.status is not None and 400 <= failure.status < 500:
            # SwanLab API 的 401/403 与其他业务拒绝：低频尝试并提示，不推断为可删除
            return "blocked"
        # 无状态码的网络异常与 5xx：可重试
        return "transient"

    # ── 调度 ──

    def schedule(self, failures: Sequence[UploadFailure], failed_rounds: int) -> RetrySchedule:
        """一轮失败 → 下次尝试间隔；多个失败取最长约束（plan §6.2）。

        随机抖动不得提前违反服务端要求：抖动只作用于退避分量，最终取各失败约束的最大值。
        解析成功的 Retry-After 回填到 failure.retry_after，供停滞诊断展示。
        空失败列表（worker 兜底吞掉原因的未预期异常）按 blocked 低频处理。
        """
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
        """指数退避：base * multiplier^(n-1)，乘抖动后限制到普通上限（plan §6.2）。"""
        exponent = max(failed_rounds, 1) - 1
        base = min(self.base_delay * (self.multiplier**exponent), self.max_delay)
        low, high = self.JITTER_RANGE
        return min(base * self._rng(low, high), self.max_delay)

    # ── Retry-After 解析 ──

    def retry_after_of(self, failure: UploadFailure) -> Optional[float]:
        """从失败的原始 ApiError 响应头解析 Retry-After，非 ApiError 失败返回 None。"""
        if not isinstance(failure.exception, ApiError):
            return None
        response = getattr(failure.exception, "response", None)
        header = getattr(response, "headers", {}).get("Retry-After")
        return self.parse_retry_after(header)

    def parse_retry_after(self, value: Optional[str]) -> Optional[float]:
        """Retry-After 头 → 相对秒数；缺失、负值、格式无效返回 None，回退普通退避。

        支持秒数与 HTTP 日期两种形式；日期用墙钟换算一次为相对时长（plan §6.2）。
        """
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
        if seconds < 0:
            return None
        return seconds


def _parse_http_date(text: str) -> Optional[float]:
    """RFC 7231 HTTP 日期字符串（如 "Wed, 21 Oct 2026 07:28:00 GMT"）→ Unix 时间戳。

    stdlib 的 email.utils.parsedate_to_datetime 是三种 HTTP 日期格式
    （IMF-fixdate / RFC 850 / asctime）的标准解析器，urllib3 自身并不解析日期形式
    的 Retry-After，故借用之。import 延迟到函数内：秒数形式（常见路径）不触发；
    且 requests 运行时本就加载了 email 包，无额外负载。解析失败返回 None；
    无时区信息时按 UTC 处理，避免用本地时区换算。
    """
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
