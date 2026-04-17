"""
@author: caddiesnew
@file: thread.py
@time: 2026/4/16
@description: 事件驱动的 Record 上传守护线程
"""

import threading
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.safe import block as safe_block

from .buffer import RecordBuffer
from .dispatch import Dispatch
from .sender import HttpRecordSender


class Transport:
    """
    事件驱动的 Record Action 线程。

    用 Condition 驱动上传事件：
    - put() 时 notify() 立刻唤醒线程，延迟接近 0
    - wait(timeout=batch_interval) 保证最大攒批间隔
    - 清空 buffer 后在锁外委托 Dispatch，不阻塞生产者
    """

    BATCH_INTERVAL: float = 1.0
    FINISH_JOIN_TIMEOUT: int = 30
    THREAD_NAME: str = "SwanLab·Transport"

    def __init__(
        self,
        batch_interval: Optional[float] = None,
        upload_callback: Optional[Callable[[int], None]] = None,
        sender: Optional[HttpRecordSender] = None,
        auto_start: bool = True,
    ):
        self._batch_interval = self.BATCH_INTERVAL if batch_interval is None else batch_interval
        self._upload_callback = upload_callback

        self._cond = threading.Condition()
        self._buf = RecordBuffer()
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._sender_closed = False

        # Transport 持有 sender，负责创建/注入/关闭
        self._sender = sender if sender is not None else HttpRecordSender()
        self._dispatcher = Dispatch(
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
        """追加 records 到 buffer 并唤醒线程。"""
        if self._finished:
            console.warning("Transport has already been finished.")
            return
        if not records:
            return
        with self._cond:
            if self._buf.extend(records) > 0:
                self._cond.notify()

    def finish(self) -> None:
        """
        通知线程停止，等待最终排空，由线程自行关闭 sender。

        设计要点：
        - _loop 的 finally 块负责 _close_sender()，保证线程退出前 sender 可用。
        - finish() 仅设置 _finished flag 并 join，不在 join 后关闭 sender，
          避免线程仍在 dispatch 时 sender 被 use-after-close。
        - join timeout 设为 30s 以覆盖弱网下多 chunk 重试场景；
          超时后线程仍为 daemon 线程会随进程退出，不会泄漏。
        """
        if self._finished:
            return
        with self._cond:
            self._finished = True
            self._cond.notify_all()
        if self._thread is None:
            # 未启动线程时直接关闭
            self._close_sender()
            return

        self._thread.join(timeout=self.FINISH_JOIN_TIMEOUT)
        if self._thread.is_alive():
            console.warning(
                "Transport thread is still running after finish timeout; "
                "it will continue retrying as a daemon and close sender on exit."
            )
        else:
            # 线程已退出，确保 sender 已关闭（_loop finally 已调用，此处为兜底）
            self._close_sender()

    def _close_sender(self) -> None:
        if self._sender_closed:
            return
        self._sender.close()
        self._sender_closed = True

    # ── 线程主循环 ──

    def _loop(self) -> None:
        pending: List[Record] = []
        network_warning_emitted = False
        try:
            while True:
                with self._cond:
                    while not pending and not self._buf and not self._finished:
                        self._cond.wait(timeout=self._batch_interval)

                    if not pending:
                        if self._buf:
                            pending = self._buf.drain()
                        elif self._finished:
                            return
                        else:
                            continue

                # 锁外委托给 Dispatch，通过返回值判断成功/失败。
                # 失败时保留 pending，避免在 drain/prepend 之间来回复制同一批 records。
                is_success, retry_records = False, pending
                with safe_block(message="Transport dispatch error"):
                    is_success, retry_records = self._dispatcher(pending)

                if is_success:
                    pending = []
                    network_warning_emitted = False
                    with self._cond:
                        if self._buf:
                            pending = self._buf.drain()
                        elif self._finished:
                            return
                else:
                    pending = retry_records
                    if not network_warning_emitted:
                        console.warning("Upload failed, network seems unavailable.")
                        network_warning_emitted = True
                    with self._cond:
                        # 失败后退避等待；notify()/notify_all() 可提前唤醒，但不会丢掉 pending。
                        self._cond.wait(timeout=self._batch_interval)
        finally:
            self._close_sender()


__all__ = [
    "Transport",
]
