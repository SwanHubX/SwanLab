"""
@author: caddiesnew
@file: utils.py
@time: 2026/3/21
@description: 上传线程通用工具
"""

import time
from queue import Queue
from typing import List


class RecordQueue:
    """线程安全的 Record 队列，支持读写权限控制。"""

    MsgType = bytes

    def __init__(self, queue: Queue, readable: bool = True, writable: bool = True):
        self._q = queue
        self._readable = readable
        self._writable = writable

    def put(self, msg: MsgType) -> None:
        if not self._writable:
            raise RuntimeError("Queue is not writable")
        self._q.put(msg)

    def put_all(self, msgs: List[MsgType]) -> None:
        if not self._writable:
            raise RuntimeError("Queue is not writable")
        for msg in msgs:
            self._q.put(msg)

    def get_all(self) -> List[MsgType]:
        if not self._readable:
            raise RuntimeError("Queue is not readable")
        msgs: List[bytes] = []
        while not self._q.empty():
            msgs.append(self._q.get())
        return msgs


class TimerFlag:
    """任务时间标识，用于判断是否到达上传间隔。"""

    def __init__(self):
        self.flag = time.time()
        self._running = True

    def can_run(self, interval: float, cancel: bool) -> bool:
        if cancel:
            return False
        if time.time() - self.flag > interval:
            self.flag = time.time()
            return True
        return False

    @property
    def running(self) -> bool:
        return self._running

    def cancel(self) -> None:
        self._running = False


__all__ = [
    "RecordQueue",
    "TimerFlag",
]
