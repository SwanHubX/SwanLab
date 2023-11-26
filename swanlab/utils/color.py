#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 19:26:22
@File: swanlab/utils/color.py
@IDE: vscode
@Description:
    终端命令行高亮类，暂时不用第三方库，使用原生方式实现，不太能兼容低版本linux
"""


# 显示格式：\033[显示方式;前景色;背景色m
# --------------------------------------------------
# 显示方式           说明
#   0             终端默认设置
#   1             高亮显示
#   4             使用下划线
#   5             闪烁
#   7             反白显示
#   8             不可见
#   22            非粗体
#   24            非下划线
#   25            非闪烁
#
# 前景色            背景色          颜色
#  30                40            黑色
#  31                41            红色
#  32                42            绿色
#  33                43            黄色
#  34                44            蓝色
#  35                45            紫红色
#  36                46            青蓝色
#  37                47            白色
# ---------------------------------------------------

colors = {
    "RED": "\033[31m",  # 红色
    "GREEN": "\033[32m",  # 绿色
    "YELLOW": "\033[33m",  # 黄色
    "BLUE": "\033[34m",  # 蓝色
    "FUCHSIA": "\033[35m",  # 紫红色
    "CYAN": "\033[36m",  # 青蓝色
    "WHITE": "\033[37m",  # 白色
}


class Colored(object):
    # ：no color
    RESET = "\033[0m"  # 终端默认颜色

    def __color_str(self, color: str, s: str) -> str:
        """私有方法，用于显示颜色字符串

        Parameters
        ----------
        color : str
            需要显示的颜色
        s : str
            字符串

        Returns
        -------
        str
            颜色字符串
        """
        return "{}{}{}".format(colors[color], s, self.RESET)

    def red(self, s: str) -> str:
        """显示红色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            红色字符串
        """
        return self.__color_str("RED", s)

    def green(self, s: str) -> str:
        """显示绿色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            绿色字符串
        """
        return self.__color_str("GREEN", s)

    def yellow(self, s: str) -> str:
        """显示黄色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            黄色字符串
        """
        return self.__color_str("YELLOW", s)

    def blue(self, s: str) -> str:
        """显示蓝色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            蓝色字符串
        """
        return self.__color_str("BLUE", s)

    def fuchsia(self, s: str) -> str:
        """显示紫红色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            紫红色字符串
        """
        return self.__color_str("FUCHSIA", s)

    def cyan(self, s: str) -> str:
        """显示青蓝色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            青蓝色字符串
        """
        return self.__color_str("CYAN", s)

    def white(self, s: str) -> str:
        """显示白色字符串

        Parameters
        ----------
        s : str
            字符串

        Returns
        -------
        str
            白色字符串
        """
        return self.__color_str("WHITE", s)


color = Colored()


if __name__ == "main":
    print(color.red("I am red!"))
    print(color.green("I am green!"))
    print(color.yellow("I am yellow!"))
    print(color.blue("I am blue!"))
    print(color.fuchsia("I am fuchsia!"))
    print(color.cyan("I am cyan!"))
    print(color.white("I am white!"))
