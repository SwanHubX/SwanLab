"""
@author: caddiesnew
@file: __init__.py
@time: 2026-03-18 00:03:03
@description: SwanLab SDK 轻量级定时任务模块

设计约束：
1. 使用守护线程执行周期任务，避免因忘记清理导致主进程无法退出
2. cancel() 为协作式停止，只会阻止下一轮调度，不会强杀正在执行的任务
3. 支持固定间隔和动态间隔策略；动态策略接收“已完成执行次数”作为输入
4. start()/run() 可重复调用，但运行中再次调用只会告警，不会启动重复线程
"""

import threading
from typing import Callable, Optional, Union

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.safe import safe_block

__all__ = ["Timer"]

Interval = Union[int, float]
IntervalStrategy = Callable[[int], Union[int, float]]
IntervalValue = Union[Interval, IntervalStrategy]


class Timer:
    """
    轻量级周期任务调度器。

    :param task: 需要周期执行的无参函数
    :param interval: 固定间隔（秒）或动态间隔策略。动态策略入参为当前已完成执行次数。
    :param immediate: 是否在启动后先立即执行一次任务，再进入下一轮等待。
    :param name: 后台线程名称，便于调试定位。
    """

    def __init__(
        self,
        task: Callable[[], None],
        *,
        interval: IntervalValue,
        immediate: bool = False,
        name: str = "SwanLab·Timer",
    ) -> None:
        self._task = task
        self._interval = interval
        self._immediate = immediate
        self._name = name

        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._count = 0

        if not callable(interval):
            self._normalize_interval(interval)

    def start(self) -> "Timer":
        """
        启动定时器。若当前已在运行，则只告警并返回自身。
        """
        with self._lock:
            if self._thread is not None and self._thread.is_alive():
                console.warning("Timer already running")
                return self

            self._stop_event.clear()
            self._thread = threading.Thread(target=self._loop, name=self._name, daemon=True)
            self._thread.start()

        return self

    def run(self) -> "Timer":
        """
        兼容旧接口，等价于 start()。
        """
        return self.start()

    def cancel(self) -> None:
        """
        发出停止信号。
        当前正在执行的任务不会被中断，但后续轮次不会再被调度。
        """
        self._stop_event.set()

    def join(self, timeout: Optional[float] = None) -> None:
        """
        等待后台线程退出，通常配合 cancel() 使用以保证任务完整收尾。
        """
        with self._lock:
            thread = self._thread

        if thread is not None and thread.is_alive():
            thread.join(timeout)

    @property
    def execution_count(self) -> int:
        """
        已完成执行次数。
        """
        return self._count

    @property
    def is_running(self) -> bool:
        """
        当前后台线程是否仍在运行。
        """
        with self._lock:
            return self._thread is not None and self._thread.is_alive()

    def _loop(self) -> None:
        if self._immediate and not self._stop_event.is_set():
            self._execute_once()

        while not self._stop_event.is_set():
            # noinspection PyBroadException
            try:
                sleep_time = self._resolve_interval()
            except Exception:
                console.trace("Timer interval strategy error")
                self._stop_event.set()
                return

            if self._stop_event.wait(sleep_time):
                return

            self._execute_once()

    def _execute_once(self) -> None:
        with safe_block(message="Error executing task"):
            self._task()
        self._count += 1

    def _resolve_interval(self) -> float:
        interval = self._interval
        if callable(interval):
            interval = interval(self._count)
        return self._normalize_interval(interval)

    @staticmethod
    def _normalize_interval(interval: Union[int, float]) -> float:
        if isinstance(interval, bool):
            raise TypeError("Timer interval must be a positive number, not bool.")

        value = float(interval)
        if value <= 0:
            raise ValueError(f"Timer interval must be greater than 0, got {interval!r}.")

        return value
