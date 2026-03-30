"""
@author: caddiesnew
@file: uploader.py
@time: 2026/3/21
@description: 上传线程池
"""

from queue import Empty, Queue
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.collector import Collector
from swanlab.sdk.internal.pkg.timer import Timer


class Uploader:
    """
    管理上传工作线程和通信管道，定时消费 BackgroundConsumer 生产的 records
    """

    UPLOAD_INTERVAL = 1.0
    UPLOAD_THREAD_NAME = "SwanLab·CorePython.Uploader"

    def __init__(
        self,
        upload_interval: Optional[float] = None,
        upload_callback: Optional[Callable] = None,
        auto_start: bool = True,
    ):
        resolved_upload_interval = self.UPLOAD_INTERVAL if upload_interval is None else upload_interval
        self._record_queue: Queue = Queue()
        self._collector = Collector(
            upload_interval=resolved_upload_interval,
            upload_callback=upload_callback,
        )
        self._timer = Timer(
            lambda: self._collector.task(self._drain_records()),
            interval=resolved_upload_interval,
            immediate=True,
            name=self.UPLOAD_THREAD_NAME,
        )
        self._started = False
        self._finished = False

        if auto_start:
            self.start()

    def start(self) -> None:
        """启动上传线程。"""
        if self._started or self._finished:
            return
        self._timer.start()
        self._started = True

    def put(self, records: List[Record]) -> None:
        """
        主线程调用：将一批 Records 投递到队列
        """
        if self._finished:
            raise RuntimeError("Uploader has already been finished.")

        self._record_queue.put(records)

    def finish(self) -> None:
        """停止线程池，执行最终 flush。"""
        if self._finished:
            return

        self._finished = True
        self._timer.cancel()
        self._timer.join(timeout=10)
        self._collector.callback(self._drain_records())

    def _drain_records(self) -> List[Record]:
        records: List[Record] = []
        while True:
            try:
                records.extend(self._record_queue.get_nowait())
            except Empty:
                return records


__all__ = [
    "Uploader",
]
