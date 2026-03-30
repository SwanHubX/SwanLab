"""
@author: caddiesnew
@file: log_collector.py
@time: 2026/3/21
@description: 上传线程聚合器
"""

import time
from typing import Callable, List, Optional

from swanlab.sdk.internal.pkg import console

from ..upload import upload_records
from .helper import RecordQueue


class UploadCollector:
    """
    日志聚合器，负责从队列中收集 Records 并批量上传。
    """

    def __init__(
        self,
        upload_interval: float = 5.0,
        upload_callback: Optional[Callable] = None,
    ):
        self.container: List[RecordQueue.MsgType] = []
        self._lock = False
        self._upload_interval = upload_interval
        self._upload_callback = upload_callback
        self._last_upload_at = time.time()

    def upload(self) -> None:
        """
        核心上传逻辑：直接批量上送 protobuf Record。
        失败的记录保留在 container 中等待下次重试。
        """
        pending = list(self.container)
        if len(pending) == 0:
            return
        upload_records(pending, upload_callback=self._upload_callback)
        self.container.clear()

    def task(self, queue: RecordQueue) -> None:
        """定时任务入口，由线程循环调用。"""
        if self._lock:
            return
        self.container.extend(queue.get_all())
        if len(self.container) == 0:
            return

        now = time.time()
        if now - self._last_upload_at <= self._upload_interval:
            return

        self._last_upload_at = now
        self._lock = True
        try:
            self.upload()
        except Exception as exc:
            console.error(f"upload error: {exc}")
        finally:
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
