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
"""
临时文件夹，每一个测试函数在执行前都会清空这个文件夹
"""

SWANLOG_FOLDER = os.path.join(TEMP_PATH, "swanlog")
"""
默认情况下，swanlog保存的文件夹
"""

SWANLAB_FOLDER = os.path.join(TEMP_PATH, ".swanlab")
"""
默认情况下，系统信息保存的文件夹
"""

# ---------------------------------- 测试用变量 ----------------------------------

TEST_CLOUD_SKIP = os.getenv("TEST_CLOUD_SKIP") is not None
"""
是否跳过云测试
"""

TEST_CLOUD_KEY = os.getenv("TEST_CLOUD_KEY")
"""
云测试的key
"""

# ---------------------------------- 导出 ----------------------------------
__all__ = [
    "TEMP_PATH",
    "nanoid",
    "TEST_CLOUD_SKIP",
    "SWANLOG_FOLDER",
    "SWANLAB_FOLDER",
    "TEST_CLOUD_KEY",
]
