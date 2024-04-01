#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/1 16:43
@File: utils.py
@IDE: pycharm
@Description:
    日志队列
"""
from queue import Queue
from typing import List, Tuple, Callable
import time

q = Queue()
"""
线程通信管道，被LogQueue类使用
"""


class LogQueue:
    """
    线程安全的日志通信队列，用于线程之间的通信
    限制了队列的可读可写性
    """

    def __init__(self, readable: bool = True, writable: bool = True):
        self.q = q
        """
        线程通信管道
        """
        self.readable = readable
        self.writable = writable

    def put(self, msg: str):
        """
        向管道中写入日志信息
        :param msg: 日志信息
        """
        if not self.writable:
            raise Exception("The queue is not writable")
        self.q.put(msg)

    def get(self):
        """
        从管道中读取日志信息
        :return: 日志信息
        """
        if not self.readable:
            raise Exception("The queue is not readable")
        return self.q.get()

    def get_all(self):
        """
        从管道中读取所有的日志信息
        :return: 日志信息
        """
        if not self.readable:
            raise Exception("The queue is not readable")
        msgs = []
        while not self.q.empty():
            msgs.append(self.q.get())
        return msgs


class TimerFlag:
    """
    任务时间标识，标识上一次任务距现在的时间，单位为秒
    可以用于判断任务是否需要执行，并且重置它
    """

    def __init__(self):
        self.flag = time.time()
        self.__running = True

    def can_run(self, interval: float, cancel: bool) -> bool:
        """
        判断任务是否可以执行
        :param interval: 任务执行间隔时间
        :param cancel: 是否强制不执行
        :return: 是否可以执行
        """
        if cancel:
            return False
        if time.time() - self.flag > interval:
            self.flag = time.time()
            return True
        return False

    @property
    def running(self):
        return self.__running

    def cancel(self):
        """
        取消任务，运行此函数后，会在下一次事件循环中退出，这代表着这个子线程也将退出
        """
        self.__running = False


class ThreadUtil:
    """
    每个线程都会传入此类的实例，统一管理线程在任务中可访问的资源
    """

    def __init__(self, queue: LogQueue, callbacks: List[Tuple[Callable, Tuple]], name: str):
        """
        初始化线程工具类
        :param queue: 线程安全的队列，用于所有线程与主要线程的通信
        :param callbacks: 回调函数列表，于让当前线程注册线程正常退出时的回调函数任务
        :param name: 线程名称
        """
        self.__queue = queue
        """
        线程安全的队列
        """
        self.__callbacks = callbacks
        self.__timer = TimerFlag()
        # 此线程是否已经注册了回调函数
        self.__set = False
        self.__name = name

    @property
    def name(self):
        return self.__name

    @property
    def queue(self) -> LogQueue:
        return self.__queue

    @property
    def timer(self):
        return self.__timer

    def register_callback(self, callback: Callable, args: Tuple):
        """
        注册回调函数
        :param callback: 回调函数
        :param args: 回调函数参数
        """
        self.__callbacks.append((callback, args))
        self.__set = True
