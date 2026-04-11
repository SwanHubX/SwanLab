"""
@author: cunyue
@file: async_task.py
@time: 2026/4/11 14:39
@description: 异步任务管理器，专为 swanlab.async_log 设计

线程安全说明：
  submit() 和 shutdown() 在 Run._api_lock 下调用，互斥。
  唯一在 api-lock 外运行的是 Future 完成回调（_on_done），回调中会从 _handles 移除已完成的 handle，
  因此 _handles 的增删需要 _handles_lock 保护。

回调时序说明：
  CPython 的 Future 实现中，set_result/set_exception 先唤醒等待者（notify_all），
  再执行回调列表（_invoke_callbacks）。因此 f.result() 可能在回调执行之前返回。
  这意味着 shutdown 仅等 f.result() 是不够的——on_success（即 Run.log()）
  可能还没执行，finish() 就先发了信号关闭 consumer，导致数据丢失。

  解决方案：每个任务关联一个 threading.Event，在回调末尾 set，shutdown 同时
  等待 future 完成和回调完成。回调完成后自动从 _handles 中移除 handle，
  保证 _handles 不会无限增长。
"""

import asyncio
import multiprocessing
import threading
import time
import traceback
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from typing import Any, Callable, Optional

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run import AsyncLogType


class _TaskHandle:
    """每个 async_log 任务关联的句柄，用于追踪回调完成状态。"""

    __slots__ = ("future", "callback_done")

    def __init__(self, future: Future):
        self.future = future
        # 回调完成后 set，shutdown 同时等待此事件
        self.callback_done = threading.Event()


class AsyncTaskManager:
    """异步任务管理器，将 func 包装为对应模式的异步任务，完成后通过回调写入日志。

    线程池、进程池均为懒初始化，首次 submit 使用对应 mode 时才创建。
    asyncio 模式复用当前运行中的事件循环，不自行创建。
    """

    def __init__(self, max_thread_workers: int = 4, max_process_workers: int = 2):
        self._max_thread_workers = max_thread_workers
        self._max_process_workers = max_process_workers

        self._thread_pool: Optional[ThreadPoolExecutor] = None
        self._spawn_pool: Optional[ProcessPoolExecutor] = None
        self._fork_pool: Optional[ProcessPoolExecutor] = None
        self._handles: list[_TaskHandle] = []
        # 保护 _handles 的增删：submit（api-lock 内）和 _on_done（api-lock 外）并发访问
        self._handles_lock = threading.Lock()

    # ----------------------------------
    # 懒初始化
    # ----------------------------------

    def _ensure_thread_pool(self) -> ThreadPoolExecutor:
        if self._thread_pool is None:
            self._thread_pool = ThreadPoolExecutor(max_workers=self._max_thread_workers)
        return self._thread_pool

    def _ensure_spawn_pool(self) -> ProcessPoolExecutor:
        if self._spawn_pool is None:
            ctx = multiprocessing.get_context("spawn")
            self._spawn_pool = ProcessPoolExecutor(max_workers=self._max_process_workers, mp_context=ctx)
        return self._spawn_pool

    def _ensure_fork_pool(self) -> ProcessPoolExecutor:
        if self._fork_pool is None:
            ctx = multiprocessing.get_context("fork")
            self._fork_pool = ProcessPoolExecutor(max_workers=self._max_process_workers, mp_context=ctx)
        return self._fork_pool

    # ----------------------------------
    # 公开接口
    # ----------------------------------

    def submit(
        self,
        func: Callable,
        args: tuple = (),
        kwargs: Optional[dict] = None,
        step: Optional[int] = None,
        mode: AsyncLogType = "threading",
        on_success: Optional[Callable[[Any, Optional[int]], None]] = None,
        on_error: Optional[Callable[[str], None]] = None,
    ) -> Future:
        """提交异步任务，返回 Future。

        :param func: 要执行的函数（asyncio 模式下须为 ``async def``）
        :param args: 位置参数
        :param kwargs: 关键字参数
        :param step: 日志 step，透传给 on_success
        :param mode: 执行模式
        :param on_success: 成功回调 ``on_success(result, step)``
        :param on_error: 失败回调 ``on_error(traceback_str)``
        :return: concurrent.futures.Future
        """
        if kwargs is None:
            kwargs = {}

        if mode == "asyncio":
            future = self._submit_asyncio(func, args, kwargs)
        elif mode == "threading":
            future = self._ensure_thread_pool().submit(func, *args, **kwargs)
        elif mode == "spawn":
            future = self._ensure_spawn_pool().submit(func, *args, **kwargs)
        elif mode == "fork":
            future = self._ensure_fork_pool().submit(func, *args, **kwargs)
        else:
            raise ValueError(f"Invalid mode: {mode!r}, expected one of 'asyncio', 'threading', 'spawn', 'fork'")

        handle = _TaskHandle(future)
        with self._handles_lock:
            self._handles.append(handle)
        future.add_done_callback(lambda f: self._on_done(f, step, on_success, on_error, handle))
        return future

    def shutdown(self, timeout: Optional[float] = None) -> None:
        """等待所有已提交任务完成并关闭资源。finish 时调用。

        等价于 wait_all + 释放线程池/进程池。调用后不可再 submit。
        timeout 为总等待时间，而非单任务等待时间。
        """
        # 1. 快照当前 handles，避免遍历时被回调修改列表
        with self._handles_lock:
            handles = list(self._handles)

        # 2. 等待所有任务执行完毕且回调已触发
        deadline = None if timeout is None else time.monotonic() + timeout
        for h in handles:
            # 2.1 检测是否超时
            remaining = None if deadline is None else max(0.0, deadline - time.monotonic())
            if remaining is not None and remaining <= 0:
                pending = sum(1 for h2 in handles if not h2.callback_done.is_set())
                if pending > 0:
                    console.warning(f"async_log: timeout reached, {pending} task(s) still pending")
                break
            # 2.2 等待任务完成
            # noinspection PyBroadException
            try:
                h.future.result(timeout=remaining)
            except Exception:
                # 任务异常已通过 on_error 回调处理，此处仅确保等待不中断
                pass
            remaining = None if deadline is None else max(0.0, deadline - time.monotonic())
            if remaining is not None and remaining <= 0:
                break
            # 等待回调完成（on_success/on_error 执行完毕）
            h.callback_done.wait(timeout=remaining)

        # 3. 释放池资源，虽然我们已经等待所有任务完成，但是更安全的做法是等待池中所有工作线程空闲后再释放
        if self._thread_pool is not None:
            self._thread_pool.shutdown(wait=True)
        if self._spawn_pool is not None:
            self._spawn_pool.shutdown(wait=True)
        if self._fork_pool is not None:
            self._fork_pool.shutdown(wait=True)

    # ----------------------------------
    # asyncio 桥接
    # ----------------------------------

    @staticmethod
    def _submit_asyncio(func: Callable, args: tuple, kwargs: dict) -> Future:
        """将协程函数提交到当前事件循环，桥接为 concurrent.futures.Future。"""
        loop = asyncio.get_running_loop()
        bridge = Future()

        async def _run():
            try:
                result = await func(*args, **kwargs)
                bridge.set_result(result)
            except Exception as e:
                bridge.set_exception(e)

        loop.create_task(_run())
        return bridge

    # ----------------------------------
    # 回调（在 api-lock 之外运行）
    # ----------------------------------

    def _on_done(
        self,
        future: Future,
        step: Optional[int],
        on_success: Optional[Callable[[Any, Optional[int]], None]],
        on_error: Optional[Callable[[str], None]],
        handle: _TaskHandle,
    ) -> None:
        # noinspection PyBroadException
        try:
            result = future.result()
            if on_success is not None:
                on_success(result, step)
        except Exception:
            if on_error is not None:
                on_error(traceback.format_exc())
        finally:
            # 从 _handles 中移除已完成的 handle，防止无限增长
            with self._handles_lock:
                self._handles.remove(handle)
            handle.callback_done.set()
