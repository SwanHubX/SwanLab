"""
@author: caddiesnew
@file: collector.py
@time: 2026/3/21
@description: 上传线程聚合器
"""

import time
from typing import List

from swanlab.sdk.internal.pkg import console

from ..upload import RecordTransport, upload_records
from .utils import RecordQueue, TimerFlag


class UploadCollector:
    """
    日志聚合器，负责从队列中收集 Records 并批量上传。
    """

    def __init__(
        self,
        transport: RecordTransport,
        upload_interval: float = 5.0,
    ):
        self._transport = transport
        self.container: List[RecordQueue.MsgType] = []
        self._lock = False
        self._upload_interval = upload_interval

    def upload(self) -> None:
        """
        核心上传逻辑：直接批量上送 protobuf Record。
        失败的记录保留在 container 中等待下次重试。
        """
        pending = list(self.container)
        if len(pending) == 0:
            return
        upload_records(pending, transport=self._transport)
        self.container.clear()

    def task(self, queue: RecordQueue, timer: TimerFlag) -> None:
        """定时任务入口，由线程循环调用。"""
        if self._lock:
            return
        self.container.extend(queue.get_all())
        if timer.can_run(self._upload_interval, len(self.container) == 0):
            self._lock = True
            try:
                self.upload()
            except Exception as exc:
                console.error(f"upload error: {exc}")
            self._lock = False

    def callback(self, queue: RecordQueue) -> None:
        """结束回调，在主线程中执行最终 flush。"""
        while self._lock:
            time.sleep(0.1)
        self.container.extend(queue.get_all())
        self.upload()


__all__ = [
    "UploadCollector",
]
