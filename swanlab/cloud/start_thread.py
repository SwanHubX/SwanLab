#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 15:58
@File: start_thread.py.py
@IDE: pycharm
@Description:
    生成线程池，生成通信管道
"""
import threading
from typing import List, Tuple, Callable, Dict
import asyncio
from .utils import ThreadUtil
from .utils import LogQueue, TimerFlag
from ._log_collector import LogCollector


class ThreadPool:
    """
    线程池类，负责管理线程和为线程生成通信管道
    只生成一个管道，用于其他线程与主要线程的通信
    每一个线程都将是一个死循环，
    """

    SLEEP_TIME = 1
    """
    默认的任务休眠时间，单位秒
    """
    MAIN_THREAD_NAME = "MsgUploader"

    def __init__(self, main_sleep_time: float = None):
        self.thread_pool = {}
        # 日志聚合器
        self.collector = LogCollector()
        # 主要线程的名称，一个线程池只能有一个主要线程
        self.main_thread = None
        # timer集合
        self.thread_timer: Dict[str, TimerFlag] = {}
        self.__callbacks: List[Tuple[Callable, Tuple]] = []
        # 生成主要线程，主要线程是日志聚合器，负责收集所有线程向主线程发送的日志信息
        self.main_thread = self.create_thread(target=self.collector.task,
                                              args=(),
                                              name=self.MAIN_THREAD_NAME,
                                              sleep_time=main_sleep_time)

    def create_thread(self,
                      target: Callable,
                      args: Tuple = (),
                      name: str = None,
                      sleep_time: float = None
                      ) -> threading.Thread:
        """
        创建一个线程
        :param target: 线程任务，约定传入的线程任务的第一个参数为 ThreadUtil 的实例
        :param args: 任务参数
        :param name: 线程名称
        :param sleep_time: 任务休眠时间
        :return: 线程对象
        """
        if name is None:
            name = f"Thread-{len(self.thread_pool)}"
        if name in self.thread_pool:
            raise Exception(f"Thread name {name} already exists")
        if sleep_time is None:
            sleep_time = self.SLEEP_TIME
        if name == self.MAIN_THREAD_NAME:
            q = LogQueue(readable=True, writable=False)
        else:
            q = LogQueue(readable=False, writable=True)
        thread_util = ThreadUtil(q, self.__callbacks, name)
        thread = threading.Thread(target=self._create_loop,
                                  args=(sleep_time, target, (thread_util, *args)),
                                  name=name)
        self.thread_pool[name] = thread
        thread.start()
        return thread

    @property
    def sub_threads(self):
        ts = []
        for name, thread in self.thread_pool.items():
            if name == self.MAIN_THREAD_NAME:
                continue
            ts.append([name, thread])
        return ts

    def _create_loop(self,
                     sleep_time: float,
                     task: Callable,
                     args: Tuple[ThreadUtil, ...]
                     ) -> asyncio.AbstractEventLoop:
        """
        创建一个事件循环，循环执行传入线程池的任务
        :param sleep_time: 任务休眠时间
        :param task: 任务
        :param args: 任务参数
        :return: 新的事件循环函数，用于启动新的线程
        """
        timer: TimerFlag = args[0].timer
        self.thread_timer[threading.current_thread().name] = timer
        # 设置事件循环为当前线程的事件循环
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        # 新的执行函数，执行任务后等待sleep_time时间后再重新执行
        async def new_task():
            while True:
                # 如果task是同步函数，直接调用，否则使用await调用
                if asyncio.iscoroutinefunction(task):
                    await task(*args)
                else:
                    task(*args)
                await asyncio.sleep(sleep_time)
                # 如果定时器停止，则退出循环
                if not timer.running:
                    return

        loop.run_until_complete(new_task())

        return loop

    def finish(self):
        """
        [在主线程中] 结束线程池中的所有线程，并执行所有线程的结束任务
        """
        print("线程池准备停止")
        # 第一步停止所有非主要线程
        for name, _ in self.sub_threads:
            self.thread_timer[name].cancel()
        for _, thread in self.sub_threads:
            thread.join()
        print("非主要线程结束")
        # 停止主要线程的任务
        self.thread_timer[self.MAIN_THREAD_NAME].cancel()
        self.main_thread.join()
        print("线程池结束")
