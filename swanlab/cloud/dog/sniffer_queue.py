#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/5 21:23
@File: sniffer_queue.py
@IDE: pycharm
@Description:
    嗅探器队列，负责收集所有嗅探线程注册的日志信息
"""
from typing import List
from queue import Queue
from ..utils import LogQueue


class SnifferQueue(LogQueue):

    def __init__(self, queue: Queue, readable: bool = True, writable: bool = True):
        super().__init__(queue, readable, writable)

    def put(self, msg: LogQueue.MsgType):
        """
        向管道中写入日志信息，日志信息必须是函数，聚合器会依次执行他们
        :param msg: 日志信息
        """
        super().put(msg)

    def get(self) -> LogQueue.MsgType:
        """
        从管道中读取日志信息
        :return: 日志信息
        """
        return super().get()

    def get_all(self) -> List[LogQueue.MsgType]:
        """
        从管道中读取所有的日志信息
        :return: 日志信息
        """
        return super().get_all()

    def put_all(self, msgs: List[LogQueue.MsgType]):
        """
        向管道中写入所有日志信息
        :param msgs: 日志信息
        """
        super().put_all(msgs)
