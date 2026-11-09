"""finish 阶段 deliver / confirm 的有限重试：只重试 retryable 失败。"""

from __future__ import annotations

import time
from typing import Callable, Optional, Protocol, TypeVar

from swanlab.sdk.protocol import CoreProtocol, CoreProtocolError

__all__ = ["FINISH_RETRY_ATTEMPTS", "FINISH_RETRY_INTERVAL", "retry_finish_call"]

FINISH_RETRY_ATTEMPTS = 3
FINISH_RETRY_INTERVAL = 1.0


class _FinishResponse(Protocol):
    success: bool
    message: str


_R = TypeVar("_R", bound=_FinishResponse)


def retry_finish_call(
    core: CoreProtocol,
    call: Callable[[], _R],
    *,
    attempts: int = FINISH_RETRY_ATTEMPTS,
    interval: float = FINISH_RETRY_INTERVAL,
    sleep: Optional[Callable[[float], None]] = None,
) -> _R:
    """执行一次 finish 调用，按需重试；返回最后一次响应。"""
    if attempts < 1:
        raise ValueError("attempts must be >= 1")
    do_sleep = time.sleep if sleep is None else sleep

    response = call()
    for _ in range(1, attempts):
        if response.success or not _is_retryable(core):
            return response
        do_sleep(interval)
        response = call()
    return response


def _is_retryable(core: CoreProtocol) -> bool:
    error = core.last_error
    return isinstance(error, CoreProtocolError) and error.retryable
