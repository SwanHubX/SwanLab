#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/28 16:33
@File: utils.py
@IDE: pycharm
@Description:
    一些工具函数
"""
from nanoid import generate
from .config import KEY


def get_password(prompt: str):
    # 如果是第一次登录，使用错误的key，会提示重新输入
    if "Paste" in prompt:
        return generate()
    else:
        return KEY
