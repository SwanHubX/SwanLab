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

from .dispatch import Dispatch
from .sender import HttpRecordSender

# 连续上传回滚超限后放弃 buffer 退出
_MAX_CONSECUTIVE_FAILURES = 3


class Transport:
    """
    事件驱动的 Record Action 线程。

    用 Condition 驱动上传事件：
    - put() 时 notify() 立刻唤醒线程，延迟接近 0
    - wait(timeout=batch_interval) 保证最大攒批间隔
    - 清空 buffer 后在锁外委托 Dispatch，不阻塞生产者
    """

    BATCH_INTERVAL = 1.0
    THREAD_NAME = "SwanLab·Transport"

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
        self._buffer: List[Record] = []
        self._buffer_num_index: set[int] = set()
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None
        self._sender_closed = False

        # Transport 持有 sender，负责创建/注入/关闭
        self._sender = sender if sender is not None else HttpRecordSender()
        self._dispatcher = Dispatch(
            cond=self._cond,
            buffer=self._buffer,
            buffer_num_index=self._buffer_num_index,
            upload_callback=self._upload_callback,
            sender=self._sender,
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
            raise RuntimeError("Transport has already been finished.")
        if not records:
            return
        with self._cond:
            filtered = [record for record in records if self._should_enqueue(record)]
            if filtered:
                self._buffer.extend(filtered)
                self._cond.notify()

    def _should_enqueue(self, record: Record) -> bool:
        """根据 record num 去重"""
        if record.num in self._buffer_num_index:
            return False
        self._buffer_num_index.add(record.num)
        return True

    def finish(self) -> None:
        """通知线程停止，等待最终排空，关闭 sender。"""
        if self._finished:
            return
        with self._cond:
            self._finished = True
            self._cond.notify_all()
        if self._thread is None:
            self._close_sender()
            return

        self._thread.join(timeout=5)
        if self._thread.is_alive():
            console.warning("Transport thread is still running after finish timeout; sender will close on thread exit.")

    def _close_sender(self) -> None:
        if self._sender_closed:
            return
        self._sender.close()
        self._sender_closed = True

    # ── 线程主循环 ──

    def _loop(self) -> None:
        consecutive_failures = 0
        try:
            while True:
                with self._cond:
                    while not self._buffer and not self._finished:
                        self._cond.wait(timeout=self._batch_interval)

                    if not self._buffer and self._finished:
                        return

                    if self._finished and consecutive_failures >= _MAX_CONSECUTIVE_FAILURES:
                        # 连续回滚超限，放弃剩余 records 退出
                        # 数据已在 DataStoreWriter 本地持久化，不影响完整性
                        console.warning(
                            f"Dropping {len(self._buffer)} records after {_MAX_CONSECUTIVE_FAILURES} consecutive upload failures"
                        )
                        self._buffer.clear()
                        self._buffer_num_index.clear()
                        return

                    pending = self._buffer[:]
                    self._buffer.clear()
                    self._buffer_num_index.clear()

                # 锁外委托给 Dispatch，通过返回值判断成功/失败
                success = self._dispatcher(pending)

                if success:
                    consecutive_failures = 0
                else:
                    consecutive_failures += 1
                    with self._cond:
                        # 失败后进入退避，但 finish() 可通过 notify_all() 立刻唤醒
                        self._cond.wait_for(lambda: self._finished, timeout=self._batch_interval)
        finally:
            self._close_sender()


__all__ = [
    "Transport",
]
