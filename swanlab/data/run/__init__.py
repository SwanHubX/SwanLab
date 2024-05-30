#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:56:48
@File: swanlab/data/run/__init__.py
@IDE: vscode
@Description:
    在此处导出SwanLabRun类，一次实验运行应该只有一个SwanLabRun实例
"""
from .main import SwanLabRun, get_run, get_config, SwanLabRunState


def register(*args, **kwargs) -> SwanLabRun:
    """注册并实例化SwanLabRun类"""
    run = SwanLabRun(*args, **kwargs)
    return run
