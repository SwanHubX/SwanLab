"""
@author: caddiesnew
@file: thread.py
@time: 2026/4/16
@description: 事件驱动的 Record 上传守护线程
"""

import threading
import time
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, safe

from .buffer import RecordBuffer
from .dispatch import Dispatch
from .sender import HttpRecordSender


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
        sender: HttpRecordSender,
        ctx: RunContext,
        upload_callback: Optional[Callable[[int], None]] = None,
        auto_start: bool = True,
    ):
        core_settings = ctx.config.settings.core
        self._batch_interval = core_settings.record.batch_interval
        self._upload_callback = upload_callback

        self._cond = threading.Condition()
        self._buf = RecordBuffer()
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._throttle = UploadWarningThrottle()

        # Transport 持有 sender，负责创建/注入/关闭
        self._sender = sender
        self._dispatcher = Dispatch(
            batch_size=core_settings.record.batch_size,
            sender=self._sender,
            upload_callback=self._upload_callback,
        )

        if auto_start:
            self.start()

    def start(self) -> None:
        """启动守护线程。"""
        if self._started or self._finished:
            return
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
            self._buf.extend(records)

    def finish(self) -> None:
        """
        通知线程停止，等待最终排空，由线程自行关闭 sender。
        """
        if self._finished:
            return
        with self._cond:
            self._finished = True
            self._cond.notify_all()
        if self._thread is None:
            return
        console.debug("Waiting for Transport to finish...")
        self._thread.join(timeout=self.FINISH_JOIN_TIMEOUT)
        console.debug("Transport finished.")

    # ── 线程主循环 ──

    def _loop(self) -> None:
        pending: List[Record] = []

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
