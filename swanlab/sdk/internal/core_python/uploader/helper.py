"""
@author: caddiesnew
@file: helper.py
@time: 2026/3/21
@description: 上传线程通用工具
"""

from queue import Queue
from typing import List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record


class RecordQueue:
    """线程安全的 Record 队列，供 uploader 生产者/消费者共享。"""

    def __init__(self, queue: Optional[Queue[Record]] = None):
        self._q: Queue[Record] = queue or Queue()

    def put(self, msg: Record) -> None:
        self._q.put(msg)

    def put_all(self, msgs: List[Record]) -> None:
        for msg in msgs:
            self._q.put(msg)

    def get(self, block: bool = True, timeout: Optional[float] = None) -> Record:
        return self._q.get(block=block, timeout=timeout)

    def get_nowait(self) -> Record:
        return self._q.get_nowait()

    def empty(self) -> bool:
        return self._q.empty()

    def qsize(self) -> int:
        return self._q.qsize()


__all__ = [
    "RecordQueue",
]
