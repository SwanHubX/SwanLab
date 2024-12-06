#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/9/20 14:35
@File: namer.py
@IDE: pycharm
@Description:
    命名器、取色器
"""
import random
from typing import Tuple, Optional

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
        prefix = random.choice(prefix_list)
        return prefix
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
        return choice_color_light, choice_color_dark  # 返回对应索引的颜色，如果超出范围则取模
