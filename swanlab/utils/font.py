#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-11 21:41:08
@File: swanlab/utils/color.py
@IDE: vscode
@Description:
    颜色处理工具
"""


import re

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


def generate_color(number: int = 1) -> str:
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
    def clear(s: str) -> str:
        """清除终端中的颜色

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


if __name__ == "__main__":
    str = """SwanLab INFO [2023-12-20 17:35:36,552] SwanLab Experiment Dashboard ready in [1m764ms

[0m[32m			➜[0m  Local:   [1mhttp://127.0.0.1:5092[0m"""
    print(FONT.clear(str))
