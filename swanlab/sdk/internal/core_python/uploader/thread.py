"""
@author: caddiesnew
@file: thread.py
@time: 2026/3/19
@description: 上传线程，定时任务 + 缓冲区 + 线程池防止数据丢失
"""

import threading
import time
from pathlib import Path
from queue import Queue
from typing import Any, Callable, Dict, List, Optional, Tuple

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console

from .model import UploadType, classify_record
from .upload import UPLOAD_DISPATCH


class RecordQueue:
    """线程安全的 Record 队列，支持读写权限控制。"""

    MsgType = Tuple[UploadType, List[Any]]

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
        msgs = []
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


class UploadCollector:
    """
    日志聚合器，负责从队列中收集 Records 并按类型聚合上传。
    对应 raw_uploader 的 LogCollectorTask。
    """

    def __init__(
        self,
        upload_interval: float = 5.0,
        upload_callback: Optional[Callable] = None,
        files_dir: Optional[Path] = None,
    ):
        self.container: List[RecordQueue.MsgType] = []
        self._lock = False
        self._upload_interval = upload_interval
        self._upload_callback = upload_callback
        self._files_dir = files_dir

    def upload(self) -> None:
        """
        核心上传逻辑：按类型聚合 → 分发到对应上传函数。
        失败的任务保留在 container 中等待下次重试。
        """
        upload_tasks: Dict[UploadType, list] = {t: [] for t in UploadType}
        for msg_type, msg_data in self.container:
            upload_tasks[msg_type].extend(msg_data)

        success_types: List[UploadType] = []
        for upload_type in UploadType:
            data_list = upload_tasks[upload_type]
            if len(data_list) == 0:
                continue
            dispatch = UPLOAD_DISPATCH.get(upload_type.value)
            if dispatch is None:
                continue
            try:
                kwargs = {}
                if dispatch["has_callback"] and self._upload_callback:
                    kwargs["upload_callback"] = self._upload_callback
                if upload_type is UploadType.FILE:
                    kwargs["files_dir"] = self._files_dir
                dispatch["upload"](data_list, **kwargs)
                success_types.append(upload_type)
            except Exception as e:
                console.error(f"{upload_type.name} upload error: {e}")

        # 失败的任务保留在 container 中
        self.container = [
            (t, upload_tasks[t]) for t in UploadType if t not in success_types and len(upload_tasks[t]) > 0
        ]

    def task(self, queue: RecordQueue, timer: TimerFlag) -> None:
        """定时任务入口，由线程循环调用。"""
        if self._lock:
            return
        self.container.extend(queue.get_all())
        if timer.can_run(self._upload_interval, len(self.container) == 0):
            self._lock = True
            try:
                self.upload()
            except Exception as e:
                console.error(f"upload error: {e}")
            self._lock = False

    def callback(self, queue: RecordQueue) -> None:
        """结束回调，在主线程中执行最终 flush。"""
        while self._lock:
            time.sleep(0.1)
        self.container.extend(queue.get_all())
        self.upload()


class ThreadPool:
    """
    上传线程池，管理上传线程和通信管道。
    保留线程池设计防止单线程数据丢失。
    """

    SLEEP_TIME = 1.0
    UPLOAD_THREAD_NAME = "SwanLab·Uploader"

    def __init__(
        self,
        upload_interval: float = 5.0,
        upload_callback: Optional[Callable] = None,
        files_dir: Optional[Path] = None,
    ):
        self._queue: Queue = Queue()
        self._collector = UploadCollector(
            upload_interval=upload_interval,
            upload_callback=upload_callback,
            files_dir=files_dir,
        )
        self._timer = TimerFlag()
        self._thread: Optional[threading.Thread] = None
        # 对外暴露只写队列
        self.queue = RecordQueue(queue=self._queue, readable=False, writable=True)

    def start(self) -> None:
        """启动上传线程。"""
        reader = RecordQueue(queue=self._queue, readable=True, writable=False)

        def loop():
            while self._timer.running:
                self._collector.task(reader, self._timer)
                time.sleep(self.SLEEP_TIME)

        self._thread = threading.Thread(target=loop, daemon=True, name=self.UPLOAD_THREAD_NAME)
        self._thread.start()

    def put(self, records: List[Record]) -> None:
        """
        主线程调用：将 Records 分类后投递到队列。
        不可上传的 Record（run/finish）会被跳过。
        """
        for record in records:
            upload_type = classify_record(record)
            if upload_type is None:
                continue
            # 将 record 序列化为 bytes 投递，由上传函数处理
            self.queue.put((upload_type, [record.SerializeToString()]))

    def finish(self) -> None:
        """停止线程池，执行最终 flush。"""
        self._timer.cancel()
        if self._thread is not None:
            self._thread.join(timeout=10)
        reader = RecordQueue(queue=self._queue, readable=True, writable=False)
        self._collector.callback(reader)


__all__ = [
    "ThreadPool",
    "RecordQueue",
    "TimerFlag",
    "UploadCollector",
]
