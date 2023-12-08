#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 19:26:22
@File: swanlab/utils/color.py
@IDE: vscode
@Description:
    颜色处理工具
"""
import random


def generate_color() -> str:
    """生成十六进制颜色字符串

    Returns
    -------
    str
        颜色字符串,以#开头的十六进制字符串,如#FFFFFF
        字符串字母大写
    """
    # 生成 RGB 随机值
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)

    # 将 RGB 转换为十六进制，并格式化为 "#RRGGBB" 的形式
    hex_color = "#{:02X}{:02X}{:02X}".format(r, g, b)

    return hex_color


if __name__ == "__main__":
    print(generate_color())
