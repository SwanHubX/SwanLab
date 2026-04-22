"""
@author: cunyue
@file: __init__.py
@time: 2026/4/22 13:56
@description: 执行器工具
"""

import threading
from concurrent.futures import Future, ThreadPoolExecutor
from typing import Any, Callable


def _is_shutting_down() -> bool:
    """检测 Python 是否正在关闭。

    threading 模块在 atexit 中设置 _SHUTTING_DOWN = True，
    此后 ThreadPoolExecutor 无法启动新线程，submit 后死锁。
    该属性为私有，但无公开 API 替代。
    """
    return getattr(threading, "_SHUTTING_DOWN", False)


class SafeThreadPoolExecutor(ThreadPoolExecutor):
    """在 Python 退出阶段仍能安全使用的线程池。

    - run(): 正常时异步提交（fire-and-forget），shutting down 时同步执行，保证任务不丢
    - submit(): shutting down 时抛 RuntimeError 而非静默死锁
    - shutdown(): 默认 wait=True，确保已提交任务执行完毕
    """

    def submit(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        """重写：shutting down 时同步执行并返回已完成 Future，而非死锁。"""
        if _is_shutting_down():
            future: Any = Future()
            try:
                future.set_result(fn(*args, **kwargs))
            except BaseException as exc:
                future.set_exception(exc)
            return future
        return super().submit(fn, *args, **kwargs)

    def run(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        """正常时异步提交，shutting down 时同步执行。

        调用方无需关心当前 Python 状态，任务一定会在某处执行。
        """
        if _is_shutting_down():
            fn(*args, **kwargs)
            return
        super().submit(fn, *args, **kwargs)

    def shutdown(self, wait: bool = True, *, cancel_futures: bool = False) -> None:
        """默认 wait=True，确保已提交任务执行完毕。"""
        super().shutdown(wait=wait, cancel_futures=cancel_futures)
