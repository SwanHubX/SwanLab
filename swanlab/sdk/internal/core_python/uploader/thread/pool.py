"""
@author: caddiesnew
@file: pool.py
@time: 2026/3/21
@description: 上传线程池
"""

import threading
import time
from queue import Queue
from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record

from ..upload import RecordTransport
from .collector import UploadCollector
from .utils import RecordQueue, TimerFlag


class ThreadPool:
    """
    上传线程池，管理上传线程和通信管道。
    保留线程池设计防止单线程数据丢失。
    """

    SLEEP_TIME = 1.0
    UPLOAD_THREAD_NAME = "SwanLab·Uploader"

    def __init__(
        self,
        transport: RecordTransport,
        upload_interval: float = 5.0,
        auto_start: bool = True,
    ):
        self._transport = transport
        self._queue: Queue = Queue()
        self._collector = UploadCollector(
            transport=transport,
            upload_interval=upload_interval,
        )
        self._timer = TimerFlag()
        self._thread: Optional[threading.Thread] = None
        self._started = False
        self._finished = False
        self.queue = RecordQueue(queue=self._queue, readable=False, writable=True)

        if auto_start:
            self.start()

    def start(self) -> None:
        """启动上传线程。"""
        if self._started or self._finished:
            return

        reader = RecordQueue(queue=self._queue, readable=True, writable=False)

        def loop() -> None:
            while self._timer.running:
                self._collector.task(reader, self._timer)
                time.sleep(self.SLEEP_TIME)

        self._thread = threading.Thread(target=loop, daemon=True, name=self.UPLOAD_THREAD_NAME)
        self._thread.start()
        self._started = True

    def put(self, records: List[Record]) -> None:
        """
        主线程调用：将 Records 序列化后投递到队列。
        """
        if self._finished:
            raise RuntimeError("ThreadPool has already been finished.")

        for record in records:
            self.queue.put(record.SerializeToString())

    def finish(self) -> None:
        """停止线程池，执行最终 flush。"""
        if self._finished:
            return

        self._finished = True
        self._timer.cancel()
        try:
            if self._thread is not None:
                self._thread.join(timeout=10)
            reader = RecordQueue(queue=self._queue, readable=True, writable=False)
            self._collector.callback(reader)
        finally:
            self._transport.close()


__all__ = [
    "ThreadPool",
]
