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

SWANLAB_LOG_DIR = os.path.join(TEMP_PATH, "swanlog")
"""
测试时swanlog文件夹存放的位置
"""

PACKAGE_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "package.mock.json")

# 注入环境变量
os.environ["SWANLAB_DEV"] = "TRUE"
os.environ["SWANLAB_PACKAGE_PATH"] = PACKAGE_PATH
os.environ["SWANLAB_LOG_DIR"] = SWANLAB_LOG_DIR

__all__ = ["TEMP_PATH", "SWANLAB_LOG_DIR", "nanoid"]
