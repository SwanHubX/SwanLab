#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/21 17:01
@File: __init__.py
@IDE: pycharm
@Description:
    tutils模块的初始化文件
"""

from .config import *
from .utils import *


def open_dev_mode() -> str:
    """
    开启开发模式，此时会返回开发环境的api-key并且创建测试目录
    在上层config部分已经执行了环境变量注入
    :return: api-key
    """
    return KEY
