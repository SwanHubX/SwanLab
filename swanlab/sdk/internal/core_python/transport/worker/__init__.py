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
from .utils import RetryPolicy, RetrySchedule


class Transport:
    """
    定时攒批的 Record Action 线程。

    用 Condition 驱动上传事件：
    - put() 仅写入 buffer，不唤醒线程
    - wait(timeout=batch_interval) 周期性唤醒线程，攒批后统一 drain
    - finish() 时 notify_all() 立刻唤醒线程排空剩余 buffer
    - 清空 buffer 后在锁外委托 Dispatch，不阻塞生产者

    失败轮的重试由 RetryPolicy 调度：指数退避 / Retry-After / 低频受阻间隔（plan §6）。
    等待可被唤醒（新入队、finish），但每次唤醒重新检查截止时间，不绕过退避。

    finish() 的返回值区分"排空成功"与"线程退出"：只有线程正常退出且缓冲与未确认
    任务均已清空（_drained）才返回 True；worker 崩溃或等待超时返回 False，
    未确认记录与最近一轮失败原因仍由实例持有（pending_records() / last_failures()
    线程退出后可查询）。
    """

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
        # 未确认任务：从取任务快照到结果应用期间始终由实例持有，worker 退出后仍可查询
        self._pending: List[Record] = []
        self._last_failures: List[UploadFailure] = []
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._throttle = UploadWarningThrottle()
        # 重试调度：失败轮截止时间（单调时钟）与连续失败计数（整轮成功后清零）
        self._retry = retry if retry is not None else RetryPolicy()
        self._retry_at: Optional[float] = None
        self._failed_rounds = 0
        # 排空标志与 worker 崩溃原因：finish() 返回值的依据（plan §8.1/§8.4）
        self._drained = False
        self._worker_error: Optional[BaseException] = None

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
        # Transport 持有实际使用的 sender，负责在 worker 退出时关闭其自持资源
        self._sender = sender
        self._dispatcher = Dispatch(
            batch_size=self._batch,
            sender=sender,
            tracker=self._tracker,
        )
        assert self._dispatcher is not None, "Dispatcher must be initialized when transport starting"
        self._thread = threading.Thread(target=self._loop, name=self.THREAD_NAME, daemon=True)
        self._thread.start()
        self._started = True

    def put(self, records: List[Record]) -> bool:
        """追加 records 到 buffer，等待下次 batch_interval 唤醒时统一处理。

        返回是否接受：draining 后的入队在锁内判定拒绝并向调用方报告（plan §5），
        不再仅打印错误；空输入视为无操作（非拒绝）。
        """
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
        """通知线程停止并等待最终排空（plan §8.1）。

        只有线程在 timeout 内退出且缓冲与未确认任务均已清空才返回 True；
        等待超时、worker 崩溃退出或从未启动但有积压数据时返回 False，
        未确认记录保留在本实例中，调用方不得据此报告上传成功。
        """
        self.request_finish()
        if self._thread is None:
            # 从未启动：无数据视为已排空，缓冲中的记录无法上传
            return not self._buf and not self._pending
        if not self.join(timeout=timeout):
            return False
        return self._drained

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

        失败重试按 RetryPolicy 计算的截止时间等待（plan §6.3）：等待可被唤醒
        （finish 等控制信号），但每次唤醒重新检查截止时间，新入队、finish 或
        重复通知都不会绕过退避或 Retry-After。
        快照在执行期间保持登记在 self._pending 上，线程意外终止时数据仍可查询。
        返回 None 表示缓冲与未确认任务均已空且已 finish，应退出。
        """
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
                    # 缓冲与未确认任务均已空且已 finish：正常排空退出
                    self._drained = True
                    return
                result: Optional[UploadResult] = None
                with safe.block(message="Transport dispatch error", write_to_tty=False):
                    result = self._dispatcher(pending)
                if result is None:
                    # 未预期异常（编程错误等）：保留本轮任务，交由节流告警与低频重试
                    remaining, failures = pending, []
                else:
                    remaining, failures = result.remaining_records, result.failures
                schedule: Optional[RetrySchedule] = None
                with self._cond:
                    self._pending = remaining
                    self._last_failures = failures
                    if not remaining and not failures:
                        # 整轮成功：清零连续失败计数与重试等待（plan §6.2）
                        self._failed_rounds = 0
                        self._retry_at = None
                    else:
                        self._failed_rounds += 1
                        schedule = self._retry.schedule(failures, self._failed_rounds)
                        self._retry_at = self._retry.clock() + schedule.delay
                if schedule is None:
                    self._throttle.reset()
                else:
                    self._throttle.warn(message=self._failure_message(schedule, failures, remaining))
        except BaseException as e:
            # worker 内部未预期异常（plan §6.1）：标记内部失败并结束线程，
            # 未确认任务保留在实例上，由 finish() 向调用方报告失败
            self._worker_error = e
            console.error(f"Transport worker exited unexpectedly: {type(e).__name__}: {e}")
        finally:
            # 关闭 sender 自持资源（清理边界，吞错不掩盖主流程）
            if self._sender is not None:
                with safe.block(message="Failed to close sender resources"):
                    self._sender.close()
            # 无论 dispatch 是否出错，transport 线程退出时标记 tracker 为 FINISHED，
            # 让进度展示停止轮询；是否排空成功由 _drained 区分。
            if self._tracker is not None:
                self._tracker.set_state(CoreState.CORE_STATE_FINISHED)

    @staticmethod
    def _failure_message(schedule: RetrySchedule, failures: List[UploadFailure], remaining: List[Record]) -> str:
        """节流告警文案：只含安全摘要（来源/阶段/状态/业务码）、重试间隔与保留数量。"""
        if failures:
            prefix = f"Upload failed ({failures[0].summary})"
        else:
            prefix = "Upload did not complete"
        return f"{prefix}; retrying in {schedule.delay:.0f}s; {len(remaining)} record(s) retained in this process."


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
