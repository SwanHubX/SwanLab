"""
@author: caddiesnew
@file: utils.py
@time: 2026/10/11
@description: worker 重试策略：指数退避
"""

import math
import random
import time
from typing import Callable, Optional, Tuple


class RetryPolicy:
    """根据上传失败信息计算下次重试等待时长（纯指数退避）。"""

    BASE_DELAY: float = 1.0
    MULTIPLIER: float = 2.0
    MAX_DELAY: float = 10.0
    JITTER_RANGE: Tuple[float, float] = (0.8, 1.2)

    def __init__(
        self,
        *,
        base_delay: Optional[float] = None,
        multiplier: Optional[float] = None,
        max_delay: Optional[float] = None,
        clock: Callable[[], float] = time.monotonic,
        rng: Optional[Callable[[float, float], float]] = None,
    ) -> None:
        self.base_delay = base_delay if base_delay is not None else self.BASE_DELAY
        self.multiplier = multiplier if multiplier is not None else self.MULTIPLIER
        self.max_delay = max_delay if max_delay is not None else self.MAX_DELAY
        self.clock = clock
        self._rng = rng if rng is not None else random.uniform
        # 指数上限，防止失败轮数过大时幂运算溢出
        if self.multiplier > 1.0 and self.max_delay > self.base_delay:
            self._saturation_exponent = math.ceil(math.log(self.max_delay / self.base_delay, self.multiplier))
        else:
            self._saturation_exponent = 0

    def schedule(self, failed_rounds: int) -> float:
        """根据连续失败轮数计算指数退避间隔。"""
        return self._backoff_delay(failed_rounds)

    def _backoff_delay(self, failed_rounds: int) -> float:
        """计算带随机抖动的指数退避等待时长。"""
        exponent = max(failed_rounds, 1) - 1
        exponent = min(exponent, self._saturation_exponent)
        base = min(self.base_delay * (self.multiplier**exponent), self.max_delay)
        low, high = self.JITTER_RANGE
        return min(base * self._rng(low, high), self.max_delay)


__all__ = [
    "RetryPolicy",
]
