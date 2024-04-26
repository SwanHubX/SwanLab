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
import shutil
import os


def clear():
    """
    清空临时文件夹
    """
    if os.path.exists(TEMP_PATH):
        shutil.rmtree(TEMP_PATH)
    os.mkdir(TEMP_PATH)
    os.mkdir(SWANLAB_LOG_DIR)


def init_db():
    """
    初始化数据库
    """
    from swanlab.db import connect, Project
    clear()
    connect(autocreate=True)
    Project.init(name="pytest-swanlab", description="测试swanlab")
