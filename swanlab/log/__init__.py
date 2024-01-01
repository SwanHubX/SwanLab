#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 17:33:56
@File: swanlab/log/__init__.py
@IDE: vscode
@Description:
    日志记录模块
"""
from typing import Optional
from .log import SwanLog

swanlog: Optional["SwanLog"] = None


def register(output_path: str, use_console_log: bool) -> SwanLog:
    """注册日志模块

    Parameters
    ----------
    output_path : str
        输出路径
    use_console_log : bool
        是否记录控制台打印信息

    Returns
    -------
    SwanLog
        初始化后的swanlog对象
    """
    global swanlog
    swanlog = SwanLog()
    swanlog.init(path=output_path, isTrain=use_console_log)
    return swanlog
