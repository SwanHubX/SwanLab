"""
@author: cunyue
@file: __init__.py
@time: 2026/3/9 19:09
@description: SwanLab 实验辅助函数，同时作为
"""

import colorsys
import random
import string
from typing import Literal, Optional, Union

__all__ = ["generate_color", "generate_id", "generate_name"]


def generate_id(length: int = 8) -> str:
    """
    Generate a unique ID for a run.

    :param length: The length of the ID.

    :return: A unique ID.
    """
    if length <= 0 or length > 64:
        raise ValueError("Length must be between 1 and 64.")
    characters = string.ascii_lowercase + string.digits
    return "".join(random.choices(characters, k=length))


PRESET_COLORS = [
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


def generate_color(slug: Optional[Union[Literal["beauty"], int]] = None) -> str:
    """
    Generate a specific or random color in hexadecimal format.

    :param slug: The slug for the color determination.
                 - If None, returns a random color from the preset list.
                 - If "beauty", generates a visually appealing random color.
                 - If an integer, returns a color from the preset list based on modulo.

    :return: A color in hexadecimal format (e.g., "#FF5733").
    """
    # 逻辑 1: 如果传入是 None，随机取一个内部列表中的值
    if slug is None:
        return random.choice(PRESET_COLORS)

    # 逻辑 2: 如果设置为 "beauty"，根据漂亮颜色算法返回一个颜色
    if slug == "beauty":
        # 漂亮颜色算法：在 HSV 空间中固定适中的饱和度(S)和较高的明度(V)，随机生成色相(H)
        h = random.random()  # 色相：0.0 到 1.0 的随机值（涵盖所有色调）
        s = random.uniform(0.4, 0.6)  # 饱和度：0.4-0.6 之间，色彩柔和不刺眼
        v = random.uniform(0.85, 1.0)  # 明度：0.85-1.0 之间，保持颜色明亮清透

        # 将 HSV 转换为 RGB (结果为 0.0 到 1.0 的浮点数)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)

        # 转换为 16 进制字符串并格式化
        return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"

    # 逻辑 3: 如果传入是 int，根据预设的内部颜色列表取模
    if isinstance(slug, int):
        return PRESET_COLORS[slug % len(PRESET_COLORS)]

    # 兜底机制（防范类型注解未生效的情况）
    return "#000000"


# 纯粹的动物名词列表
PRESET_ANIMALS = [
    "swan",
    "rat",
    "ox",
    "tiger",
    "rabbit",
    "dragon",
    "snake",
    "horse",
    "goat",
    "monkey",
    "rooster",
    "dog",
    "pig",
    "cat",
    "elephant",
    "penguin",
    "kangaroo",
    "panda",
    "lion",
    "zebra",
]

# "beauty" 模式下的优美形容词
BEAUTY_ADJECTIVES = [
    "elegant",
    "stellar",
    "vibrant",
    "serene",
    "cosmic",
    "lucid",
    "radiant",
    "ethereal",
    "nimble",
    "valiant",
    "gentle",
    "brilliant",
    "luminous",
    "tranquil",
    "dazzling",
]


def generate_name(slug: Optional[Union[Literal["beauty"], int]] = None) -> str:
    """
    Generate a specific or random entity name.

    :param slug: The slug for the name determination.
                 - If None, returns a random animal + 4-char random hash (e.g., "swan-a3f9").
                 - If "beauty", generates an appealing adjective-animal-number combo (e.g., "stellar-swan-42").
                 - If an integer, returns a name based on modulo + the integer itself (e.g., "dragon-128").

    :return: A generated name string.
    """
    # 逻辑 1: 如果传入是 None，随机动物 + 4位随机字符后缀
    if slug is None:
        animal = random.choice(PRESET_ANIMALS)
        random_hash = "".join(random.choices(string.ascii_lowercase + string.digits, k=4))
        return f"{animal}-{random_hash}"

    # 逻辑 2: 如果设置为 "beauty"，生成形容词+动物+数字的组合
    if slug == "beauty":
        adj = random.choice(BEAUTY_ADJECTIVES)
        animal = random.choice(PRESET_ANIMALS)
        number = random.randint(10, 99)
        return f"{adj}-{animal}-{number}"

    # 逻辑 3: 如果传入是 int，根据预设列表取模，并将 ID 拼在尾部
    if isinstance(slug, int):
        animal = PRESET_ANIMALS[slug % len(PRESET_ANIMALS)]
        return f"{animal}-{slug}"

    # 兜底机制
    fallback_hash = "".join(random.choices(string.ascii_lowercase + string.digits, k=4))
    return f"unknown-{fallback_hash}"
