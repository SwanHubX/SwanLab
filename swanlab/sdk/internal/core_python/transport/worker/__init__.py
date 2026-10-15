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
from .utils import RetryPolicy


class Transport:
    """周期攒批并上传 Record 的后台守护线程。"""

    THREAD_NAME: str = "SwanLab·Transport"
    FINISH_JOIN_TIMEOUT: int = 30

    def __init__(
        self,
        ctx: CoreContext,
        tracker: Optional[UploadTracker] = None,
        sender: Optional[HttpRecordSender] = None,
        auto_start: bool = True,
        retry: Optional[RetryPolicy] = None,
    ):
        self._ctx = ctx
        self._batch = ctx.config.record_batch
        self._batch_interval = ctx.config.record_interval
        self._tracker = tracker
        self._sender = sender

        self._cond = threading.Condition()
        self._buf = RecordBuffer()
        # 未确认 records
        self._pending: List[Record] = []
        self._last_failures: List[UploadFailure] = []
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._throttle = UploadWarningThrottle()
        # 重试调度
        self._retry = retry if retry is not None else RetryPolicy()
        self._retry_at: Optional[float] = None
        self._failed_rounds = 0
        # 排空标志与 worker 崩溃原因
        self._drained = False
        self._worker_error: Optional[BaseException] = None

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
        self._sender = sender
        self._dispatcher = Dispatch(
            batch_size=self._batch,
            sender=sender,
            tracker=self._tracker,
        )
        self._thread = threading.Thread(target=self._loop, name=self.THREAD_NAME, daemon=True)
        self._thread.start()
        self._started = True

    def put(self, records: List[Record]) -> bool:
        """追加 records 到缓冲区，返回是否接受入队。"""
        if not records:
            return True
        accepted = 0
        with self._cond:
            if self._finished:
                rejected = True
            else:
                rejected = False
                accepted = self._buf.extend(records)
        if rejected:
            console.error("Transport has already been finished; records were rejected.")
            return False
        if self._tracker is not None:
            self._tracker.add_total_records(accepted)
        return True

    def request_finish(self) -> None:
        """通知线程停止并排空缓冲区。"""
        with self._cond:
            if self._finished:
                return
            self._finished = True
            self._cond.notify_all()

    def join(self, timeout: Optional[float] = FINISH_JOIN_TIMEOUT) -> bool:
        """等待线程退出，返回是否在 timeout 内退出。"""
        if self._thread is None:
            return True
        console.debug("Waiting for Transport to finish...")
        self._thread.join(timeout=timeout)
        console.debug("Transport finished.")
        return not self._thread.is_alive()

    def is_alive(self) -> bool:
        """后台线程是否运行中。"""
        return self._thread is not None and self._thread.is_alive()

    def finish(self, timeout: Optional[float] = FINISH_JOIN_TIMEOUT) -> bool:
        """通知线程停止并等待缓冲区与未确认记录排空。"""
        self.request_finish()
        if self._thread is None:
            return not self._buf and not self._pending
        if not self.join(timeout=timeout):
            return False
        return self._drained

    def pending_records(self) -> List[Record]:
        """返回未确认记录与缓冲区记录快照。"""
        with self._cond:
            return list(self._pending) + self._buf.snapshot()

    def last_failures(self) -> List[UploadFailure]:
        """返回最近一轮上传失败列表快照。"""
        with self._cond:
            return list(self._last_failures)

    def finish_summary(self) -> str:
        """生成排空失败的单行诊断摘要。"""
        reasons: List[str] = []
        if self._worker_error is not None:
            reasons.append(f"worker crashed: {type(self._worker_error).__name__}")
        failures = self.last_failures()
        if failures:
            reasons.append(f"last failure: {failures[0].summary}")
        if not reasons:
            reasons.append("worker exited unexpectedly")
        retained = len(self.pending_records())
        return f"Upload did not finish: {retained} record(s) retained ({'; '.join(reasons)})."

    def _take_snapshot(self) -> Optional[List[Record]]:
        """获取下一批待上传记录，失败重试按截止时间等待。"""
        with self._cond:
            if self._pending:
                deadline = self._retry_at
                if deadline is not None:
                    remaining = deadline - self._retry.clock()
                    while remaining > 0:
                        self._cond.wait(timeout=remaining)
                        deadline = self._retry_at
                        if deadline is None:
                            break
                        remaining = deadline - self._retry.clock()
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
                    self._drained = True
                    return
                result: Optional[UploadResult] = None
                with safe.block(message="Transport dispatch error", write_to_tty=False):
                    result = self._dispatcher(pending)
                if result is None:
                    remaining, failures = pending, []
                else:
                    remaining, failures = result.remaining_records, result.failures
                retry_delay: Optional[float] = None
                with self._cond:
                    self._pending = remaining
                    self._last_failures = failures
                    if not remaining:
                        self._failed_rounds = 0
                        self._retry_at = None
                    else:
                        self._failed_rounds += 1
                        retry_delay = self._retry.schedule(self._failed_rounds)
                        self._retry_at = self._retry.clock() + retry_delay
                if retry_delay is None:
                    self._throttle.reset()
                else:
                    self._throttle.warn()
        except BaseException as e:
            self._worker_error = e
            console.error(f"Transport worker exited unexpectedly: {type(e).__name__}: {e}")
        finally:
            if self._sender is not None:
                with safe.block(message="Failed to close sender resources"):
                    self._sender.close()
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)


class UploadWarningThrottle:
    """控制上传失败警告输出频率，恢复时提示。"""

    INTERVAL: float = 30.0

    def __init__(self) -> None:
        # 初始值设为 -INTERVAL，使首次调用触发告警
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
        self._last_warn_at = -self.INTERVAL
        self._in_failure = False


__all__ = [
    "Transport",
]
