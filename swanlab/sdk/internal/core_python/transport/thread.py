"""
@author: caddiesnew
@file: thread.py
@time: 2026/4/16
@description: 事件驱动的 Record 上传守护线程
"""

import threading
import time
from typing import Any, List, Optional

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.pkg import console, safe

from .buffer import RecordBuffer
from .dispatch import Dispatch
from .sender import HttpRecordSender
from .tracker import UploadTracker


class Transport:
    """
    定时攒批的 Record Action 线程。

    用 Condition 驱动上传事件：
    - put() 仅写入 buffer，不唤醒线程
    - wait(timeout=batch_interval) 周期性唤醒线程，攒批后统一 drain
    - finish() 时 notify_all() 立刻唤醒线程排空剩余 buffer
    - 清空 buffer 后在锁外委托 Dispatch，不阻塞生产者
    """

    THREAD_NAME: str = "SwanLab·Transport"
    FINISH_JOIN_TIMEOUT: int = 30

    def __init__(
        self,
        ctx: CoreContext,
        tracker: Optional[UploadTracker] = None,
        sender: Optional[Any] = None,
        auto_start: bool = True,
        track_record_totals: bool = True,
    ):
        self._ctx = ctx
        self._batch = ctx.config.record_batch
        self._batch_interval = ctx.config.record_interval
        self._tracker = tracker
        self._sender = sender
        self._track_record_totals = track_record_totals

        self._cond = threading.Condition()
        self._buf = RecordBuffer()
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._throttle = UploadWarningThrottle()

        # Transport 持有 sender，负责创建/注入/关闭
        self._dispatcher: Optional[Dispatch] = None
        if auto_start:
            self.start()

    def start(self) -> None:
        """启动守护线程。"""
        if self._started or self._finished:
            return
        sender = self._sender or HttpRecordSender(ctx=self._ctx)
        if self._tracker is not None:
            sender.set_tracker(self._tracker)
        self._dispatcher = Dispatch(
            batch_size=self._batch,
            sender=sender,
        )
        assert self._dispatcher is not None, "Dispatcher must be initialized when transport starting"
        self._thread = threading.Thread(target=self._loop, name=self.THREAD_NAME, daemon=True)
        self._thread.start()
        self._started = True

    def put(self, records: List[Record]) -> None:
        """追加 records 到 buffer，等待下次 batch_interval 唤醒时统一处理。"""
        if self._finished:
            console.error("Transport has already been finished.")
            return
        if not records:
            return
        with self._cond:
            accepted = self._buf.extend(records)
        if self._tracker is not None and self._track_record_totals:
            self._tracker.add_total_records(accepted)

    def request_finish(self) -> None:
        """通知线程停止并排空 buffer，不等待线程退出。"""
        if self._finished:
            return
        with self._cond:
            self._finished = True
            self._cond.notify_all()

    def join(self, timeout: Optional[float] = FINISH_JOIN_TIMEOUT) -> bool:
        """等待线程退出，返回是否在 timeout 内正常退出。"""
        if self._thread is None:
            return True
        console.debug("Waiting for Transport to finish...")
        self._thread.join(timeout=timeout)
        console.debug("Transport finished.")
        return not self._thread.is_alive()

    def is_alive(self) -> bool:
        """Transport worker thread is still running."""
        return self._thread is not None and self._thread.is_alive()

    def finish(self, timeout: Optional[float] = FINISH_JOIN_TIMEOUT) -> bool:
        """
        通知线程停止，等待最终排空，由线程自行关闭 sender。
        """
        self.request_finish()
        return self.join(timeout=timeout)

    # ── 线程主循环 ──

    def _loop(self) -> None:
        pending: List[Record] = []
        assert self._dispatcher is not None, "Dispatcher must be initialized when transport starting"

        def refill_pending() -> bool:
            """尝试为 pending 补充数据。

            :returns
                True: 可以继续循环
                False: 没有更多数据且已 finish，应退出
            """
            nonlocal pending

            with self._cond:
                while not self._buf and not self._finished:
                    self._cond.wait(timeout=self._batch_interval)
                if self._buf:
                    pending = self._buf.drain()
                    return True
                return not self._finished

        try:
            while True:
                if not pending and not refill_pending():
                    return
                is_success, retry_records = False, pending
                with safe.block(message="Transport dispatch error"):
                    is_success, retry_records = self._dispatcher(pending)
                if is_success:
                    pending = []
                    self._throttle.reset()
                else:
                    pending = retry_records
                    self._throttle.warn()
        finally:
            # 无论 dispatch 是否出错，transport 线程退出时标记 tracker 为 FINISHED，
            # 让进度展示停止轮询。
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)


class UploadWarningThrottle:
    """控制上传失败警告的打印频率，成功后重置并提示恢复。"""

    INTERVAL: float = 30.0

    def __init__(self) -> None:
        # 用 -INTERVAL 而非 0.0，确保 time.monotonic() 返回值小于 INTERVAL 时
        # （如 CI 容器刚启动）仍能首次触发 warn()
        self._last_warn_at: float = -self.INTERVAL
        self._in_failure: bool = False

    def warn(self, message: str = "Upload failed due to network or server issues, retrying automatically.") -> None:
        self._in_failure = True
        now = time.monotonic()
        if now - self._last_warn_at >= self.INTERVAL:
            console.warning(message)
            self._last_warn_at = now

    def reset(self, message: str = "Upload recovered, resuming normally.") -> None:
        if self._in_failure:
            console.info(message)
        # 同 __init__，用 -INTERVAL 保证 reset 后下次 warn() 立刻触发
        self._last_warn_at = -self.INTERVAL
        self._in_failure = False


__all__ = [
    "Transport",
]
