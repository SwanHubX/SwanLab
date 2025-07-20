"""
@author: cunyue
@file: log.py
@time: 2025/5/15 18:35
@description: 标准输出、标准错误流拦截代理，支持外界设置/取消回调，基础作用为输出日志
"""

import re
import sys
from typing import List, Tuple, Callable

from swanlab.toolkit import SwanKitLogger, create_time, LogContent
from .counter import AtomicCounter
from .type import LogHandler, LogType, WriteHandler, LogData, ProxyType


class SwanLog(SwanKitLogger):
    """
    swanlab 日志类
    继承自 SwanKitLogger 的同时增加标准输出、标准错误留拦截代理功能
    """

    def __init__(self, name=__name__.lower(), level="info"):
        super().__init__(name=name, level=level)
        self.__original_level = level
        # 当前已经代理的输出行数
        self.__counter = AtomicCounter(0)
        # 保存原始的标准输出和标准错误流
        self.__origin_stdout_write = None
        self.__origin_stderr_write = None
        # 代理缓冲区
        self.__stdout_buffer = ""
        self.__stderr_buffer = ""
        # 上传到云端的最大长度
        self.__max_upload_len = None
        # 当前的代理类型
        self.__proxy_type = None

    @property
    def epoch(self):
        return self.__counter.value

    def __create_write_handler(self, write_type: LogType, handler: LogHandler) -> WriteHandler:
        """
        创建一个新的处理器
        """
        origin_write_handler = self.__origin_stdout_write if write_type == 'stdout' else self.__origin_stderr_write

        def get_buffer():
            return self.__stdout_buffer if write_type == 'stdout' else self.__stderr_buffer

        def set_buffer(buffer):
            if write_type == 'stdout':
                self.__stdout_buffer = buffer
            else:
                self.__stderr_buffer = buffer

        max_output_len = self.__max_upload_len

        def write_handler(message: str):
            """
            处理器函数，线程安全
            """
            try:
                origin_write_handler(message)
            except UnicodeEncodeError:
                # 遇到编码问题，直接pass，此时表现为终端不输出
                pass
            except ValueError as e:
                if "I/O operation on closed file" in str(e):
                    # 遇到文件已关闭问题，直接pass，此时表现为终端不输出
                    pass

            # 进行缓冲处理，主要目的是处理进度条输出：
            # 1. 如果 message 不包含换行符，则加入缓冲区，否则跳转步骤2
            # 2. 如果 message 包含换行符，则将换行符之前的内容和当前缓冲区合并，进入步骤3，准备上传
            # 3. 根据换行符分隔为一个个 message，解析 message，如果 message 包含\r
            messages, new_buffer = clean_control_chars(get_buffer() + message)
            set_buffer(new_buffer)
            # 4. 遍历 messages，上传到云端
            if len(messages):
                log_data = LogData(
                    type=write_type,
                    contents=[],
                )
                with self.__counter as counter:
                    for message in messages:
                        log_data['contents'].append(
                            LogContent(
                                message=message[:max_output_len],
                                create_time=create_time(),
                                epoch=counter.increment(),
                            )
                        )
                # 设置回调
                handler(log_data)

        return write_handler

    def __exec_fun_by_type(self, stdout_func: Callable, stderr_func: Callable):
        """
        根据设置的类型执行对应的函数
        :return: None
        """
        if self.__proxy_type == 'all':
            stdout_func()
            stderr_func()
        elif self.__proxy_type == 'stdout':
            stdout_func()
        elif self.__proxy_type == 'stderr':
            stderr_func()
        # 为 none 就不管了

    @property
    def proxied(self):
        """
        判断是否已经开启代理了
        :return: bool
        """
        return self.__origin_stderr_write is not None or self.__origin_stdout_write is not None

    def start_proxy(self, proxy_type: ProxyType, max_log_length: int, handler: LogHandler, epoch: int = None):
        """
        启动代理
        :param max_log_length: 一行日志的最大长度，超过这个长度的日志将被截断，-1 表示不限制
        :param proxy_type: 代理类型，支持 "stdout", "stderr", "all"
        :param handler: 代理处理函数
        :param epoch: 可选参数，设置当前起始 epoch，默认为 None
        """
        if self.proxied:
            raise RuntimeError("Std Proxy is already started")
        # 设置一些状态
        self.__max_upload_len = max_log_length
        self.__proxy_type = proxy_type
        if epoch is not None:
            self.__counter = AtomicCounter(epoch)

        # 设置代理
        def set_stdout():
            self.__stdout_buffer = ""
            self.__origin_stdout_write = sys.stdout.write
            sys.stdout.write = self.__create_write_handler('stdout', handler)

        def set_stderr():
            self.__stderr_buffer = ""
            self.__origin_stderr_write = sys.stderr.write
            sys.stderr.write = self.__create_write_handler('stderr', handler)

        self.__exec_fun_by_type(set_stdout, set_stderr)

    def stop_proxy(self):
        """
        停止代理
        """
        # 如果没有开启代理，则直接返回
        if not self.proxied:
            return

        # 清理标准输出
        def clean_stdout():
            if self.__stdout_buffer:
                sys.stdout.write(self.__stdout_buffer + '\n')
            self.__stdout_buffer = ""
            sys.stdout.write = self.__origin_stdout_write
            self.__origin_stdout_write = None

        # 清理标准错误
        def clean_stderr():
            if self.__stderr_buffer:
                sys.stderr.write(self.__stderr_buffer + '\n')
            self.__stderr_buffer = ""
            sys.stderr.write = self.__origin_stderr_write
            self.__origin_stderr_write = None

        self.__exec_fun_by_type(clean_stdout, clean_stderr)
        self.__counter = AtomicCounter(0)

    def reset(self):
        """
        重置输出流代理和日志等级
        """
        self.stop_proxy()
        self.level = self.__original_level


_ANSI_ESCAPE_RE = re.compile(r'\x1b\[[0-9;]*[a-zA-Z]')  # 匹配ANSI控制码


def clean_control_chars(text) -> Tuple[List[str], str]:
    """
    清理终端控制字符（模拟终端覆盖行为）
    特别，如果当前字符串中不包含换行符，则直接返回，不处理
    如果存在行为空字符串，则将其删除，不在结果中返回
    """
    lines = text.split('\n')
    cleaned_lines = []
    # 最后一行不处理，因为可能是未完成的行
    # 这与当前 “如果当前字符串中不包含换行符，则直接返回，不处理” 的设计是一致的，并且有更多的灵活性: 最后一行将被作为 buffer 缓冲
    for line in lines[:-1]:
        # 移除控制序列
        cleaned_line = remove_control_sequences(line)
        # 使用预编译正则处理字符串，删除ANSI控制码
        cleaned_line = _ANSI_ESCAPE_RE.sub('', cleaned_line)
        cleaned_line and cleaned_lines.append(cleaned_line)
    return cleaned_lines, lines[-1]


def remove_control_sequences(line: str) -> str:
    """
    高效移除控制序列，目前只处理 \r 和 \x1b[A：
    1. \r：回车符，表示光标回到行首
    2. \x1b[A：上移一行，通常用于终端覆盖输出
    我们取最后一个控制序列后的内容
    """
    # 查找最后一个控制序列的位置
    last_cr = line.rfind('\r')
    last_up = line.rfind('\x1b[A')
    # 确定最后出现的位置
    last_control = max(last_cr, last_up)
    if last_control == -1:
        return line  # 没有控制序列，直接返回

    # 根据控制序列类型确定截取位置
    if last_control == last_up:
        return line[last_control + 3 :]  # \x1b[A 长3字符
    else:
        return line[last_control + 1 :]  # \r 长1字符
