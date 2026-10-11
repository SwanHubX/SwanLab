"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/16
@description: 事件驱动的 Record 上传守护线程：线程、调度、入队、排空
"""

import threading
import time
from typing import List, Optional

from swanlab.deprecated.sender import DeprecatedHttpRecordSender
from swanlab.proto.swanlab.operation.v1.operation_pb2 import CoreState
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

from ..dispatch import Dispatch
from ..sender import HttpRecordSender
from ..tracker import UploadTracker
from .buffer import RecordBuffer


class Transport:
    """
    定时攒批的 Record Action 线程。

    用 Condition 驱动上传事件：
    - put() 仅写入 buffer，不唤醒线程
    - wait(timeout=batch_interval) 周期性唤醒线程，攒批后统一 drain
    - finish() 时 notify_all() 立刻唤醒线程排空剩余 buffer
    - 清空 buffer 后在锁外委托 Dispatch，不阻塞生产者

    未确认记录与最近一轮失败原因始终由实例持有（self._pending / self._last_failures），
    线程退出后仍可通过 pending_records() / last_failures() 查询。
    """

    THREAD_NAME: str = "SwanLab·Transport"
    FINISH_JOIN_TIMEOUT: int = 30

    def __init__(
        self,
        ctx: CoreContext,
        tracker: Optional[UploadTracker] = None,
        sender: Optional[HttpRecordSender] = None,
        auto_start: bool = True,
    ):
        self._ctx = ctx
        self._batch = ctx.config.record_batch
        self._batch_interval = ctx.config.record_interval
        self._tracker = tracker
        self._sender = sender

        self._cond = threading.Condition()
        self._buf = RecordBuffer()
        # 未确认任务：从取任务快照到结果应用期间始终由实例持有，worker 退出后仍可查询
        self._pending: List[Record] = []
        self._last_failures: List[UploadFailure] = []
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

        # 初始化sender，根据不同的项目版本选择不同的请求器，以兼容历史版本
        # NOTE: 相关兼容逻辑将在v0.11版本后移除
        sender = self._sender
        if sender is None:
            sender = (
                HttpRecordSender(ctx=self._ctx)
                if self._ctx._project_version == 1
                else DeprecatedHttpRecordSender(ctx=self._ctx)
            )

        if self._tracker is not None:
            sender.set_tracker(self._tracker)
        self._dispatcher = Dispatch(
            batch_size=self._batch,
            sender=sender,
            tracker=self._tracker,
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
        if self._tracker is not None:
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

    def pending_records(self) -> List[Record]:
        """当前未确认任务与缓冲中的记录快照（诊断用，线程退出后仍可查询）。"""
        with self._cond:
            return list(self._pending) + self._buf.snapshot()

    def last_failures(self) -> List[UploadFailure]:
        """最近一轮上传的失败原因快照（诊断用，线程退出后仍可查询）。"""
        with self._cond:
            return list(self._last_failures)

    # ── 线程主循环 ──

    def _take_snapshot(self) -> Optional[List[Record]]:
        """取下一批任务快照：优先重试未确认记录，其次排空缓冲。

        快照在执行期间保持登记在 self._pending 上，线程意外终止时数据仍可查询。
        返回 None 表示缓冲与未确认任务均已空且已 finish，应退出。
        """
        with self._cond:
            if self._pending:
                return self._pending
            while not self._buf and not self._finished:
                self._cond.wait(timeout=self._batch_interval)
            if self._buf:
                pending = self._buf.drain()
                self._pending = pending
                return pending
            return None

    def _loop(self) -> None:
        assert self._dispatcher is not None, "Dispatcher must be initialized when transport starting"

        try:
            while True:
                pending = self._take_snapshot()
                if pending is None:
                    return
                result: Optional[UploadResult] = None
                with safe.block(message="Transport dispatch error", write_to_tty=False):
                    result = self._dispatcher(pending)
                if result is None:
                    # 未预期异常（编程错误等）：保留本轮任务，交由节流告警
                    remaining, failures = pending, []
                else:
                    remaining, failures = result.remaining_records, result.failures
                with self._cond:
                    self._pending = remaining
                    self._last_failures = failures
                if not remaining and not failures:
                    self._throttle.reset()
                else:
                    self._throttle.warn(message=self._failure_message(failures, remaining))
        finally:
            # 无论 dispatch 是否出错，transport 线程退出时标记 tracker 为 FINISHED，
            # 让进度展示停止轮询。
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)

    @staticmethod
    def _failure_message(failures: List[UploadFailure], remaining: List[Record]) -> str:
        """节流告警文案：只含安全摘要（来源/阶段/状态/业务码）与保留数量。"""
        if failures:
            return f"Upload failed ({failures[0].summary}); {len(remaining)} record(s) retained and will be retried."
        return f"Upload did not complete; {len(remaining)} record(s) retained and will be retried."


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
