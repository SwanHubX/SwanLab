#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-11 21:41:08
@File: swanlab/utils/color.py
@IDE: vscode
@Description:
    颜色处理工具
"""

import sys
import re
from typing import Callable, Tuple
import threading
import time

light_colors = [
    "#528d59",  # 绿色
    "#587ad2",  # 蓝色
    "#c24d46",  # 红色
    "#9cbe5d",  # 青绿色
    "#6ebad3",  # 天蓝色
    "#dfb142",  # 橙色
    "#6d4ba4",  # 紫色
    "#8cc5b7",  # 淡青绿色
    "#892d58",  # 紫红色
    "#40877c",  # 深青绿色
    "#d0703c",  # 深橙色
    "#d47694",  # 粉红色
    "#e3b292",  # 淡橙色
    "#b15fbb",  # 浅紫红色
    "#905f4a",  # 棕色
    "#989fa3",  # 灰色
]

COLOR_LIST = {
    "light": light_colors,
    "dark": light_colors,
}


def generate_color(number: int = 1):
    """输入数字，在设定好顺序的颜色列表中返回十六进制颜色字符串

    Returns
    -------
    str
        颜色字符串,以#开头的十六进制字符串,如#FFFFFF
        字符串字母大写
    """

    if number % 16 == 0:
        number = 16
    else:
        number = number % 16

    return COLOR_LIST["light"][number - 1], COLOR_LIST["dark"][number - 1]


class FONT:

    @staticmethod
    def loading(
        s: str,
        func: Callable,
        args: Tuple = (),
        interval: float = 0.4,
        prefix: str = None,
        brush_length: int = 100
    ):
        """
        实现终端打印的加载效果，输入的字符串会在开头出现loading效果以等待传入的函数执行完毕

        Parameters
        ----------
        s : str
            需要打印的字符串
        func : coroutine
            执行的同步函数
        args : Tuple, optional
            传入函数的参数，默认为空
        interval : float, optional
            loading的速度，即每个字符的间隔时间，单位为秒
        prefix : str, optional
            前缀字符串，打印在loading效果之前，默认为swanlab
        brush_length : int, optional
            刷去的长度，默认为100
        """
        # FIXME 因为协程有一些适配性问题，暂时使用线程
        prefix = FONT.bold(FONT.blue("swanlab")) + ": " if prefix is None else prefix
        symbols = ["\\", "|", "/", "-"]

        running, result, error = True, None, None

        def loading():
            index = 0
            while True:
                sys.stdout.write("\r" + prefix + symbols[index % len(symbols)] + " " + s)
                sys.stdout.flush()
                index += 1
                time.sleep(interval)
                if not running:
                    break

        # 再次封装传入的func，当func执行完毕以后，将running置为False
        def task():
            nonlocal result, error, running
            try:
                result = func(*args)
            except Exception as e:
                error = e
            finally:
                running = False

        # 开启新线程
        t1 = threading.Thread(target=loading)
        t2 = threading.Thread(target=task)
        t1.start()
        t2.start()
        try:
            t2.join()
        except KeyboardInterrupt:
            running = False
            raise KeyboardInterrupt
        t1.join()
        if error is not None:
            raise error
        FONT.brush("", brush_length)
        return result

    @staticmethod
    def swanlab(s: str, color: str = "blue"):
        """用于为某一条信息添加swanlab前缀"""
        return FONT.bold(getattr(FONT, color)("swanlab")) + ": " + s

    @staticmethod
    def brush(s: str, length: int = 20) -> None:
        """
        将当前终端行刷去，替换为新的字符串

        Parameters
        ----------
        s : str
            需要刷去的字符串
        length : int, optional
            需要刷去的长度，默认为20，如果当前行的长度大于length，但又需要刷去整行，则需要传入更大的length
        """
        sys.stdout.write("\r" + " " * length + "\r" + s)
        sys.stdout.flush()

    @staticmethod
    def bold(s: str) -> str:
        """在终端中加粗字符串

        Parameters
        ----------
        s : str
            需要加粗的字符串

        Returns
        -------
        str
            加粗后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[1m{s}\033[0m"

    @staticmethod
    def default(s: str) -> str:
        """在终端中将字符串着色为默认颜色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[39m{s}\033[0m"

    @staticmethod
    def blue(s: str) -> str:
        """在终端中将字符串着色为蓝色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[34m{s}\033[0m"

    @staticmethod
    def grey(s: str) -> str:
        """在终端中将字符串着色为灰色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[90m{s}\033[0m"

    @staticmethod
    def underline(s: str) -> str:
        """在终端中将字符串着色为下划线

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[4m{s}\033[0m"

    @staticmethod
    def green(s: str) -> str:
        """在终端中将字符串着色为绿色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[32m{s}\033[0m"

    @staticmethod
    def dark_green(s: str) -> str:
        """在终端中将字符串着色为深绿色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[38;5;22m{s}\033[0m"

    @staticmethod
    def dark_gray(s: str) -> str:
        """在终端中将字符串着色为深灰色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[38;5;236m{s}\033[0m"

    @staticmethod
    def yellow(s: str) -> str:
        """在终端中将字符串着色为黄色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[33m{s}\033[0m"

    @staticmethod
    def red(s: str) -> str:
        """在终端中将字符串着色为红色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[31m{s}\033[0m"

    @staticmethod
    def magenta(s: str) -> str:
        """在终端中将字符串着色为品红色

        Parameters
        ----------
        s : str
            需要着色的字符串

        Returns
        -------
        str
            着色后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return f"\033[35m{s}\033[0m"

    @staticmethod
    def clear(s: str) -> str:
        """清除字符串中的颜色编码

        Parameters
        ----------
        s : str
            需要清除颜色的字符串

        Returns
        -------
        str
            清除颜色后的字符串
        """
        ansi_escape_pattern = re.compile(r"\033\[[0-9;]+m")
        return ansi_escape_pattern.sub("", s)
