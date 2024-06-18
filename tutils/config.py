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

__test_path = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__)
        )
    ),
    "test"
)

TEMP_PATH = os.path.join(__test_path, "temp")

SWANLAB_DIR = os.path.join(TEMP_PATH, ".swanlab")
"""
测试时.swanlab文件夹存放的位置
"""

SWANLAB_LOG_DIR = os.path.join(TEMP_PATH, "swanlog")
"""
测试时swanlog文件夹存放的位置
"""

CONFIG: dict = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.json")))
"""
开发快捷配置
"""
KEY: str = CONFIG["api-key"]
"""
测试时使用的api-key
"""

PACKAGE_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "package.mock.json")

# 注入环境变量
os.environ["SWANLAB_DEV"] = "TRUE"
os.environ["SWANLAB_PACKAGE_PATH"] = PACKAGE_PATH
os.environ["SWANLAB_LOG_DIR"] = SWANLAB_LOG_DIR
os.environ["SWANLAB_HOME"] = TEMP_PATH

__all__ = ["TEMP_PATH", "SWANLAB_LOG_DIR", "CONFIG", "nanoid", "PACKAGE_PATH", "SWANLAB_DIR", "KEY"]
