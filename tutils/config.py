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

# ---------------------------------- 路径配置 ----------------------------------

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

PACKAGE_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "package.mock.json")

SWANLOG_FOLDER = os.path.join(TEMP_PATH, "swanlog")

SWANLAB_FOLDER = os.path.join(TEMP_PATH, ".swanlab")


def reset_env():
    os.environ[SwanLabEnv.SWANLAB_PACKAGE.value] = PACKAGE_PATH
    os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = SWANLOG_FOLDER
    os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = SWANLAB_FOLDER


reset_env()

__all__ = ["TEMP_PATH", "nanoid", "KEY", "SWANLOG_FOLDER", "SWANLAB_FOLDER", "PACKAGE_PATH", "reset_env"]

# ---------------------------------- 测试用变量 ----------------------------------

SkipTest = None
"""
在某些情况下，我们需要跳过测试，这个变量用于跳过某些测试
"""
