#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 16:38:02
@File: swanlab/database/expriments_name.py
@IDE: vscode
@Description:
    随机生成实验名称，以树木的名称命名
"""
import random

# 树木名称列表
tree_names = [
    "oak",  # 橡树
    "maple",  # 枫树
    "pine",  # 松树
    "beech",  # 榉树
    "willow",  # 柳树
    "birch",  # 桦树
    "spruce",  # 云杉
    "cypress",  # 柏树
    "cedar",  # 雪松
    "hemlock",  # 铁杉
    "poplar",  # 杨树
    "redwood",  # 红杉
]

# 正向形容词列表
positive_adjectives = [
    "lush",  # 茂盛的
    "majestic",  # 壮观的
    "vibrant",  # 蓬勃的
    "verdant",  # 翠绿的
    "beautiful",  # 美丽的
    "sturdy",  # 稳固的
    "abundant",  # 富饶的
    "tall",  # 高大的
    "shady",  # 树荫浓密的
    "hardy",  # 耐寒的
]


def generate_random_tree_name(exists_name: [str] = []) -> str:
    """随机生成树木名称, 以树木的名称命名实验

    Parameters
    ----------
    exists_name : [str], optional
        已存在的树木名称, by default []

    Returns
    -------
    str
        随机生成的树木名称, 例如: lush-oak-1
    """
    # 随机选择一个树木名称和一个正向形容词
    tree_name = random.choice(tree_names)
    adjective = random.choice(positive_adjectives)

    # 拼接生成的树木名称
    generated_name = f"{adjective}-{tree_name}-{exists_name.__len__() + 1}"
    return generated_name


def make_experiment_name_unique(name: str, exists_name: [str]) -> str:
    """检查实验名称是否已经存在，如果存在，调用一个递归函数，后缀加1

    Parameters
    ----------
    name : str
        待检查的实验名称

    exists_name : [str]
        已经存在的实验名称列表

    Returns
    -------
    str
        检查后的实验名称
    """
    if name in exists_name:
        # 如果名称已经存在，递归调用自己，后缀加1，需要判断后缀是否存在以及是否可以转换为数字
        name_arr = name.split("-")
        # 如果最后一个元素是数字，那么后缀加1
        if name_arr[-1].isdigit():
            name_arr[-1] = str(int(name_arr[-1]) + 1)
        else:
            # 如果最后一个元素不是数字，那么添加一个1
            name_arr.append("1")
        # 重新拼接名称
        new_name = "-".join(name_arr)
        return make_experiment_name_unique(new_name, exists_name)
    return name


def check_experiment_name(name: str):
    """检查实验名称是否符合规范，如果不符合规范，报错

    Parameters
    ----------
    name : str
        名称
    """
    return
