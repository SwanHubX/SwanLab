#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:28
@File: config.py.py
@IDE: pycharm
@Description:
    存储一些与单元测试有关的快捷配置
"""
import os
import json
import nanoid
from swanlab.env import SwanLabEnv

__test_path = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__)
        )
    ),
    "test"
)

TEMP_PATH = os.path.join(__test_path, "temp")

_ = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.json")))

KEY: str = _["api-key"]
"""
测试时使用的api-key
"""

_ = os.path.join(os.path.abspath(os.path.dirname(__file__)), "package.mock.json")

# 注入环境变量
os.environ[SwanLabEnv.SWANLAB_PACKAGE.value] = _
os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = os.path.join(TEMP_PATH, "swanlog")
os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = os.path.join(TEMP_PATH, ".swanlab")

__all__ = ["TEMP_PATH", "nanoid", "KEY"]
