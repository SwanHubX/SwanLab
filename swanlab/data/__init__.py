#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 16:08:15
@File: swanlab/data/__init__.py
@IDE: vscode
@Description:
    在此处完成回调注册、swanlog注册，并为外界提供api，提供运行时生成的配置
"""

import atexit, sys, traceback, os
from datetime import datetime
from ..log import swanlog
from .sdk import (
    init,
    log,
    finish,
)


# def log(data: dict, step: int = None):
#     """以字典的形式记录数据，字典的key将作为列名，value将作为记录的值
#     例如:
#     ```python
#     sw.log({"loss": 0.1, "accuracy": 0.9})
#     ```
#     Parameters
#     ----------
#     data : dict
#         此处填写需要记录的数据
#     step: int
#         当前记录的步数，如果不传则默认当前步数为'已添加数据数量+1'
#     """
#     if sd is None:
#         raise RuntimeError("swanlab is not initialized")
#     if not isinstance(data, dict):
#         raise TypeError("log data must be a dict")
#     if step is not None and (not isinstance(step, int) or step < 0):
#         raise TypeError("'step' must be an integer not less than zero.")
#     for key in data:
#         # 遍历字典的key，记录到本地文件中
#         d = data[key]
#         # 检查key的类型
#         check_key_format(key)
#         # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
#         sd.add(key, d, step=step)


__all__ = ["log", "init"]
