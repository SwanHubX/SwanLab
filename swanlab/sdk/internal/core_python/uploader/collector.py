"""
@author: caddiesnew
@file: collector.py
@time: 2026/3/21
@description: 上传线程聚合器
"""

import threading
import time
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.uploader.sender import upload_records
from swanlab.sdk.internal.pkg import console


class Collector:
    """
    日志聚合器，负责从队列中收集 Records 并批量上传。
    """

    def __init__(
        self,
        upload_interval: float = 1.0,
        upload_callback: Optional[Callable] = None,
    ):
        self.container: List[Record] = []
        self._lock = threading.Lock()
        self._upload_interval = upload_interval
        self._upload_callback = upload_callback
        self._last_upload_at = time.time()

    def upload(self, pending: List[Record]) -> None:
        """
        核心上传逻辑：直接批量上送 protobuf Record。
        由于网络请求耗时，应在无锁状态下调用。
        """
        if not pending:
            return
        upload_records(pending, upload_callback=self._upload_callback)

    def task(self, records: List[Record]) -> None:
        """定时任务入口，由线程循环调用。"""
        with self._lock:
            self.container.extend(records)
            if not self.container:
                return

            now = time.time()
            if now - self._last_upload_at <= self._upload_interval:
                return

            self._last_upload_at = now
            pending = self.container
            self.container = []

        try:
            self.upload(pending)
        except Exception as exc:
            console.error(f"upload error: {exc}")
            with self._lock:
                self.container = pending + self.container

    def callback(self, records: List[Record]) -> None:
        """结束回调，在主线程中执行最终 flush。"""
        with self._lock:
            self.container.extend(records)
            pending = self.container
            self.container = []

        try:
            self.upload(pending)
        except Exception as exc:
            console.error(f"upload error: {exc}")
            with self._lock:
                self.container = pending + self.container


__all__ = [
    "Collector",
]
