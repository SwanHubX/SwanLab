"""
@author: caddiesnew
@file: thread.py
@time: 2026/4/16
@description: 事件驱动的 Record Action 线程
"""

import threading
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

from .dispatch import Dispatch


class Transport:
    """
    事件驱动的 Record Action 线程。

    用 Condition 驱动写事件：
    - put() 时 notify() 立刻唤醒线程，延迟接近 0
    - wait(timeout=batch_interval) 保证最大攒批间隔
    - swap buffer 后锁外委托 Dispatch，不阻塞生产者
    """

    BATCH_INTERVAL = 1.0
    THREAD_NAME = "SwanLab·Transport"

    def __init__(
        self,
        batch_interval: Optional[float] = None,
        upload_callback: Optional[Callable[[int], None]] = None,
        auto_start: bool = True,
    ):
        self._batch_interval = self.BATCH_INTERVAL if batch_interval is None else batch_interval
        self._upload_callback = upload_callback

        self._cond = threading.Condition()
        self._buffer: List[Record] = []
        self._finished = False
        self._started = False
        self._thread: Optional[threading.Thread] = None

        self._dispatcher = Dispatch(cond=self._cond, buffer=self._buffer, upload_callback=self._upload_callback)

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
        with self._cond:
            self._buffer.extend(records)
            self._cond.notify()

    def finish(self) -> None:
        """通知线程停止，等待最终排空。"""
        if self._finished:
            return
        with self._cond:
            self._finished = True
            self._cond.notify_all()
        if self._thread is not None:
            self._thread.join(timeout=10)

    # ── 线程主循环 ──

    def _loop(self) -> None:
        while True:
            with self._cond:
                while not self._buffer and not self._finished:
                    self._cond.wait(timeout=self._batch_interval)

                if not self._buffer and self._finished:
                    return

                pending = self._buffer
                self._buffer = []

            # 锁外委托给 Dispatch
            self._dispatcher(pending)


__all__ = [
    "Transport",
]
