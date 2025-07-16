#!/usr/bin/env python
r"""
@DATE: 2024/9/20 14:35
@File: namer.py
@IDE: pycharm
@Description:
    命名器、取色器
"""

import random
from typing import Optional, Tuple

prefix_list = [
    "swan-",  # 天鹅
    "rat-",  # 鼠
    "ox-",  # 牛
    "tiger-",  # 虎
    "rabbit-",  # 兔
    "dragon-",  # 龙
    "snake-",  # 蛇
    "horse-",  # 马
    "goat-",  # 羊
    "monkey-",  # 猴
    "rooster-",  # 鸡
    "dog-",  # 狗
    "pig-",  # 猪
    "cat-",  # 猫
    "elephant-",  # 象
    "penguin-",  # 企鹅
    "kangaroo-",  # 袋鼠
    "monkey-",  # 猴
    "panda-",  # 熊猫
    "tiger-",  # 虎
    "lion-",  # 狮
    "zebra-",  # 斑马
]


def generate_name(index: Optional[int] = None) -> str:
    """
    生成名称
    :param index: 生成名称的索引，约定为历史实验数，可以为None
    :return: 生成的名称
    """
    if index is None:
        prefix: str = random.choice(prefix_list)
        return prefix[:-1]
    else:
        prefix = prefix_list[index % len(prefix_list)]
        return prefix + str(index + 1)


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

# 黑夜模式的颜色暂时和light_colors保持一致
dark_colors = [
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


def hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
    """将十六进制颜色转换为RGB格式

    Args:
        hex_color: 十六进制颜色字符串，支持以下格式：
            - "#528d59"
            - "528d59"
            - "#fff"  (简写形式)
            - "fff"

    Returns:
        Tuple[int, int, int]: RGB格式的颜色值，范围[0,255]

    Raises:
        ValueError: 当输入格式不正确时

    Examples:
        >>> hex_to_rgb("#528d59")
        (82, 141, 89)
        >>> hex_to_rgb("fff")  # 简写形式
        (255, 255, 255)
    """
    # 去除可能存在的#号和空白
    hex_color = hex_color.strip().lstrip("#").strip()

    # 验证输入
    if not all(c in "0123456789abcdefABCDEF" for c in hex_color):
        raise ValueError(f"Invalid hex color: {hex_color}")

    # 处理简写形式 (例如 #fff)
    if len(hex_color) == 3:
        hex_color = "".join(c * 2 for c in hex_color)
    elif len(hex_color) != 6:
        raise ValueError(f"Invalid hex color length: {len(hex_color)}, should be 3 or 6 characters")

    try:
        # 将16进制转换为RGB
        r = int(hex_color[:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:], 16)

        return (r, g, b)
    except ValueError as e:
        raise ValueError(f"Invalid hex color format: {hex_color}") from e


def generate_colors(index: Optional[int] = None) -> Tuple[str, str]:
    """
    生成颜色
    :param index: 生成颜色的索引，约定为历史实验数，可以为None
    :return: 生成的颜色，(白天颜色，夜晚颜色)
    """

    if index is None:
        choice_color = random.choice(light_colors)
        return choice_color, choice_color  # 随机返回一个颜色
    else:
        choice_color_light = light_colors[index % len(light_colors)]
        choice_color_dark = dark_colors[index % len(dark_colors)]
        return (
            choice_color_light,
            choice_color_dark,
        )  # 返回对应索引的颜色，如果超出范围则取模


def generate_run_id() -> str:
    """
    生成运行ID
    :return: 生成的运行ID，21个字符的小写字母和数字
    """
    return "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=21))
