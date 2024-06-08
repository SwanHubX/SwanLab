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
from typing import Tuple, Callable, List
import time
from abc import ABC, abstractmethod
from .task_types import UploadType


class LogQueue:
    """
    日志队列，定义一些工具和可读可写性质
    """
    MsgType = Tuple[UploadType, List]
    """
    传入日志聚合器的日志信息类型，应该是一个元组，第一个元素是文件类型，第二个元素是日志信息，日志信息应该是一个列表
    """

    def __init__(self, queue: Queue, readable: bool = True, writable: bool = True):
        self.q = queue
        self.__readable = readable
        self.__writable = writable

    @property
    def readable(self):
        return self.__readable

    @property
    def writable(self):
        return self.__writable

    def put(self, msg: MsgType):
        """
        向管道中写入日志信息，日志信息必须是函数，聚合器会依次执行他们
        :param msg: 日志信息
        """
        if not self.writable:
            raise Exception("The queue is not writable")
        self.q.put(msg)

    def get(self) -> MsgType:
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

    def put_all(self, msgs: List[MsgType]):
        """
        向管道中写入所有日志信息
        :param msgs: 日志信息
        """
        if not self.writable:
            raise Exception("The queue is not writable")
        for msg in msgs:
            self.q.put(msg)


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

    def __init__(self, queue: LogQueue, name: str):
        """
        初始化线程工具类
        :param queue: 线程安全的队列，用于所有线程与主要线程的通信
        :param name: 线程名称
        """
        self.__queue = queue
        """
        线程安全的队列
        """
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

    @staticmethod
    def wrapper_callback(func: Callable, args: Tuple) -> Callable:
        """
        回调函数包装器，将回调函数和参数包装成一个函数
        :param func: 回调函数
        :param args: 回调函数参数
        """

        def wrapper():
            func(*args)

        return wrapper


class ThreadTaskABC(ABC):
    """
    线程抽象类，规定接入线程池所必须的方法和接口
    """

    @abstractmethod
    def task(self, u: ThreadUtil, **kwargs):
        """线程正常执行时的任务"""
        pass

    @abstractmethod
    def callback(self, u: ThreadUtil, *args):
        """线程执行完毕(被关闭）后的回调函数"""
        pass
