"""
@author: cunyue
@file: timer.py
@time: 2025/12/30 22:04
@description: SwanLab 封装的定时器工具，用于定时、循环执行若干函数。

采用 Daemon 线程 + Event 机制，支持手动 cancel 和 join，保证任务完整性。
"""

import threading
from typing import Callable, Union, Optional

from swanlab.log import swanlog


class Timer:
    def __init__(self, task: Callable, *, interval: Union[int, float, Callable[[int], float]], immediate: bool = False):
        """
        初始化定时器

        :param task: 需要定时执行的任务函数
        :param interval: 执行间隔（秒）。可以是固定的数字，也可以是返回数字的函数（用于动态间隔），输入为task调用、执行次数
        :param immediate 是否立即执行，默认否
        """
        self._task = task
        self._interval = interval
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._immediate = immediate

        # 互斥锁，用于标记“任务正在运行中”
        self._run_lock = threading.Lock()
        # 标记调用次数
        self._count = 0

    def run(self) -> "Timer":
        """
        启动定时器
        """
        if self._thread is not None and self._thread.is_alive():
            return swanlog.warning("Timer already running")

        # 重置停止信号，允许重启
        self._stop_event.clear()

        # daemon=True 防止忘记 cancel 导致进程挂死
        # 如果用户需要保证数据不丢失，应手动调用 cancel() + join()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        return self

    def cancel(self):
        """
        发出停止信号。
        注意：这不会强制中断正在执行的任务，而是等待当前任务执行完毕后不再进行下一次循环。
        """
        self._stop_event.set()

    def join(self, timeout=None):
        """
        等待定时器线程结束。
        通常在调用 cancel() 后调用此方法，以确保最后一次任务完整执行。
        """
        if self._thread is not None and self._thread.is_alive():
            self._thread.join(timeout)

    def _loop(self):
        """
        线程主循环
        """
        # 第一次立即执行
        if self._immediate and not self._stop_event.is_set():
            self._safe_execute()

        while not self._stop_event.is_set():
            # wait 既起到 sleep 的作用，又能响应 set() 事件
            # 如果在 sleep 期间调用了 cancel()，这里会立即唤醒并返回 True，从而跳出循环
            self._stop_event.wait(self._sleep_time)

            if not self._stop_event.is_set():
                self._safe_execute()

    def _safe_execute(self):
        """
        执行任务，并使用锁保护
        """
        with self._run_lock:
            try:
                self._task()
            except Exception as e:
                swanlog.error(f"Error executing task: {e}")
            finally:
                self._count += 1

    @property
    def _sleep_time(self) -> float:
        """
        解析间隔时间，支持动态策略
        """
        if callable(self._interval):
            return self._interval(self._count)
        return float(self._interval)
