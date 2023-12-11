#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-11 21:41:08
@File: swanlab/utils/color.py
@IDE: vscode
@Description:
    颜色处理工具
"""
import random


def rgb_to_hex(rgb_color: tuple):
    """将RGB转为十六进制颜色字符串

    Returns
    -------
    str
        颜色字符串,以#开头的十六进制字符串,如#FFFFFF
        字符串字母大写
    """
    r, g, b = rgb_color

    # 将rgb转为颜色字符串
    hex_color = "#{:02X}{:02X}{:02X}".format(r, g, b)

    return hex_color


def hex_to_rgb(hex_color: str):
    """将十六进制颜色字符串转为RGB

    Returns
    -------
    tuple
        包含RGB的元组,每个元素都是0-255之间的整数,如(255, 255, 255)
    """

    # 去除可能包含的 '#' 符号
    hex_color = hex_color.lstrip("#")

    # 将十六进制颜色代码分成红、绿和蓝部分
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    return (r, g, b)


def generate_color() -> str:
    """生成十六进制颜色字符串

    Returns
    -------
    str
        颜色字符串,以#开头的十六进制字符串,如#FFFFFF
        字符串字母大写
    """

    # 生成 RGB 随机变化值
    r_random = random.randint(0, 10)
    g_random = random.randint(0, 10)
    b_random = random.randint(0, 10)

    # 生成随机数, 用于在颜色列表中选择一个随机颜色
    random_number = random.randint(0, 15)

    color_list = [
        "#528d59",
        "#9cbe5d",
        "#dfb142",
        "#d0703c",
        "#e3b292",
        "#c24d46",
        "#892d58",
        "#d47694",
        "#8cc5b7",
        "#40877c",
        "#6ebad3",
        "#587ad2",
        "#6d4ba4",
        "#b15fbb",
        "#905f4a",
        "#989fa3",
    ]

    # 将随机选择的十六进制字符串转为RGB
    r, g, b = hex_to_rgb(color_list[random_number])

    # 在RGB通道增加随机波动
    r = min(r + r_random, 255)
    g = min(g + g_random, 255)
    b = min(b + b_random, 255)

    # 将RGB转为十六进制字符串，然后返回
    return rgb_to_hex(rgb_color=(r, g, b))


if __name__ == "__main__":
    print(generate_color())
