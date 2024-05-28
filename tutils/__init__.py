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
import shutil
import os


def clear():
    """
    清空临时文件夹, 重新创建
    """
    if os.path.exists(TEMP_PATH):
        shutil.rmtree(TEMP_PATH)
    os.mkdir(TEMP_PATH)
    os.mkdir(SWANLAB_LOG_DIR)
    os.mkdir(SWANLAB_DIR)


def init_db():
    """
    初始化数据库
    """
    from swanlab.db import connect, Project
    clear()
    connect(autocreate=True)
    Project.init(name="pytest-swanlab", description="测试swanlab")


def open_dev_mode() -> str:
    """
    开启开发模式，此时会返回开发环境的api-key并且创建测试目录
    :return: api-key
    """
    # 创建测试目录
    os.makedirs(SWANLAB_LOG_DIR, exist_ok=True)
    return CONFIG['api-key']
