"""
@author: cunyue
@file: __init__.py
@time: 2026/4/22 13:56
@description: 执行器工具
"""

import asyncio
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from concurrent.futures import thread as futures_thread
from typing import Any, Callable, Coroutine, Optional

__all__ = ["SafeThreadPoolExecutor", "EventLoopExecutor"]


def _is_shutting_down() -> bool:
    """检测 Python 是否正在关闭。

    threading 模块在 atexit 中设置 _SHUTTING_DOWN = True，
    concurrent.futures.thread 也会设置 _shutdown = True。
    此后 ThreadPoolExecutor 无法启动新线程，submit 后可能抛 RuntimeError。
    这些属性为私有，但无公开 API 替代。
    """
    return getattr(threading, "_SHUTTING_DOWN", False) or getattr(futures_thread, "_shutdown", False)


class SafeThreadPoolExecutor(ThreadPoolExecutor):
    """在 Python 退出阶段仍能安全使用的线程池。

    - run(): 正常时异步提交（fire-and-forget），shutting down 时同步执行，保证任务不丢
    - submit(): shutting down 时同步执行并返回已完成 Future
    - shutdown(): 默认 wait=True，确保已提交任务执行完毕
    """

    def submit(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        """重写：shutting down 时同步执行并返回已完成 Future，而非死锁或抛错。"""
        if _is_shutting_down():
            future: Any = Future()
            try:
                # 来自此issue: https://github.com/SwanHubX/SwanLab/issues/889，此时需要一个个发送
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


class EventLoopExecutor:
    """Run one coroutine in a dedicated background event loop."""

    def __init__(self, name: str) -> None:
        self._name = name
        self._loop: Optional[asyncio.AbstractEventLoop] = None
        self._thread: Optional[threading.Thread] = None
        self._future: Optional[Future[None]] = None

    def start(self, coro: Coroutine[Any, Any, None]) -> None:
        if self._future is not None:
            raise RuntimeError(f"Background task {self._name!r} is already running")

        loop = asyncio.new_event_loop()
        thread = threading.Thread(target=loop.run_forever, name=self._name, daemon=True)
        thread.start()

        self._loop = loop
        self._thread = thread
        self._future = asyncio.run_coroutine_threadsafe(coro, loop)

    def wait(self) -> None:
        try:
            if self._future is not None:
                self._future.result()
        finally:
            self.close()

    def close(self) -> None:
        if self._loop is not None:
            loop = self._loop
            loop.call_soon_threadsafe(loop.stop)

        if self._thread is not None:
            self._thread.join()

        if self._loop is not None:
            self._loop.close()

        self._loop = None
        self._thread = None
        self._future = None
