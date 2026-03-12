"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:29
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理
2. 运行、实验上下文管理
3. 触发回调
"""

import queue
import threading
from functools import cached_property, singledispatchmethod
from pathlib import Path
from typing import Any, Mapping, Optional, Tuple, Union, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.data.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import FinishType

from . import metrics
from .callbackers import CloudCallback, LocalCallback, OfflineCallback
from .data.transforms import Text

__all__ = ["SwanLabRun", "CloudCallback", "LocalCallback", "OfflineCallback"]

from .helper import flatten_dict

LogData = Mapping[str, Any]
"""
日志数据代表的是单次 swanlab.log({...}, step=N) 调用中的数据字典（展开后）。
"""
LogPayload = Tuple[LogData, int, Timestamp]
"""
日志负载代表的是单次 swanlab.log({...}, step=N) 调用，0代表log的日志字典（展开后），1代表step，2代表时间戳。
"""


class SwanLabRun:
    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx
        self._queue: queue.Queue[Union[LogPayload, None]] = queue.Queue(maxsize=100_000)
        # 启动后台处理线程，相关的落盘、回调触发等处理不会阻塞主线程
        self._background_logger = threading.Thread(target=self._background_log, daemon=False)
        self._background_logger.start()

    @cached_property
    def id(self) -> str:
        """
        This property returns the unique identifier of the current SwanLab run.
        """
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def run_dir(self) -> Path:
        """
        This property returns the directory where the current SwanLab run is saved.
        """
        assert self._ctx.run_dir is not None, "Run dir is not set."
        return self._ctx.run_dir

    @cached_property
    def _flush_file(self) -> Path:
        """
        指向存储文件
        """
        assert self.run_dir is not None, "Run dir is not set when attempting to access _flush_file."
        return self._ctx.backup_file

    def log(self, data: LogData, step: Optional[int] = None):
        """
        Log a row of data to the current run. Unlike `swanlab.log`, this api will be called directly based on the
        SwanRun instance, removing the initialization process. Of course, after you call the finish method,
        this method will be banned from calling.

        :param data: The data to log.
        :param step: The global step at which the data was logged. Can be None if not explicitly tracked.
        """
        # 类型检查
        if not isinstance(data, Mapping):
            console.error(f"Log data must be a dict, but got {type(data).__name__}. SwanLab will ignore it.")
            return
        # 如果没有传递step，获取全局step，并将全局step+1
        step = metrics.next_step(self._ctx, step)
        # 获取当前时间
        ts = Timestamp()
        ts.GetCurrentTime()
        # 展开字典，例如 {"a": {"b": {"c": 1}}} -> {"a/b/c": 1}
        data = flatten_dict(data)
        # 放入队列
        # 如果队列已满，阻塞等待直到有空位，否则会丢弃数据
        self._queue.put((data, step, ts), block=True)

    def _background_log(self, flush_timeout: float = 0.5, batch_size: int = 100) -> None:
        """将数据放入队列，异步处理"""
        batch_records = []
        self._registered_columns = set()
        while True:
            try:
                # 阻塞等待数据到来
                payload = self._queue.get(timeout=flush_timeout)
                # 收到 finish 信号
                if payload is None:
                    self._flush_to_disk(batch_records)
                    break
                # 解析数据
                data, step, timestamp = payload
                # 处理数据
                for key, value in data.items():
                    try:
                        # 内存态转换
                        _: LogRecord = self._dispatch_log(value, key=key, timestamp=timestamp, step=step)
                        # TODO 如果是首次写入，需要创建列

                        # TODO: 写入文件

                        # TODO: 触发回调器
                    except Exception as e:
                        console.error(f"Error in background logger: {e}")
                if len(batch_records) >= batch_size:
                    self._flush_to_disk(batch_records)
                    batch_records.clear()
            except queue.Empty:
                # 超时，把缓存的数据刷盘
                if batch_records:
                    self._flush_to_disk(batch_records)
                    batch_records.clear()
            except Exception as e:
                # 无论发生什么，后台线程不能死
                console.error(f"SwanLab background logger thread error: {e}")

    def _flush_to_disk(self, records: list):
        """
        将日志写入本地文件，落盘成功后，会通知回调器，进入下一步处理
        :param records: 待落盘的日志记录
        :return: None
        """
        if not records:
            return
        try:
            # TODO: 批量追加到本地文件
            # TODO: 落盘成功后，才通知回调器
            ...
        except Exception as e:
            console.error(f"SwanLab failed to write disk: {e}")

    @singledispatchmethod
    def _dispatch_log(self, value: Any, key: str, timestamp: Timestamp, step: Optional[int]) -> LogRecord:
        """默认的回退处理逻辑
        此部分用于处理基础标量数据类型
        """
        ...

    @_dispatch_log.register
    def _(self, value: Text, key: str, timestamp: Timestamp, step: Optional[int]) -> LogRecord:
        """Log一个Text"""
        ...

    def log_text(self, key: str, data: Union[str, Text], caption: Optional[str] = None, step: Optional[int] = None):
        """
        Log a text record to the current run.
        :param key: The key of the record.
        :param data: The text content or a Text object.
        :param caption: An optional caption for the text.
        :param step: The global step at which the data was logged. Can be None if not explicitly tracked.
        """
        ...

    def define_media(self):
        """
        Define a media record to the current run.
        """
        ...

    def define_scalar(self):
        """
        Define a scalar record to the current run.
        """
        ...

    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """
        Finish the current run and close the current experiment
        Normally, swanlab will run this function automatically,
        but you can also execute it manually and mark the experiment as 'completed'.
        Once the experiment is marked as 'completed', no more data can be logged to the experiment by 'swanlab.log'.
        After calling this function, you can re-run `swanlab.init` to start a new experiment.

        :param state: The state of the experiment, it can be 'success', 'crashed' or 'aborted'.
        :param error: The error message when the experiment is marked as 'crashed'. If not 'crashed', it should be None.
        """
        # 类型检查
        state = state.lower()  # type: ignore
        if state not in get_args(FinishType):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return
        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"
        # 清理副作用
        # 1. 清理log队列
        self._queue.put(None)
        self._background_logger.join()
