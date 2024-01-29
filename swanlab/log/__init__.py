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

swanlog: Optional["SwanLog"] = SwanLog("SwanLab")


def register(
    output_path: str = None,
    console_path: str = None,
    log_level: str = None,
    enable_logging: bool = False,
) -> SwanLog:
    """注册日志模块

    Parameters
    ----------
    output_path : str
        输出路径
    use_console_log : str
        是否记录控制台打印信息
    log_level : str
        日志等级
    enable_logging : bool
        是否启用日志

    Returns
    -------
    SwanLog
        初始化后的swanlog对象
    """
    swanlog.init(path=output_path, console_path=console_path, level=log_level, enable_logging=enable_logging)
    return swanlog
