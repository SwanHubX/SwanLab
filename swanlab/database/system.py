#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:49:47
@File: swanlab\database\system.py
@IDE: vscode
@Description:
    采集系统数据，包括内存、CPU、GPU、硬盘、网络等
"""
import platform
import socket


def get_system_info():
    """获取系统信息"""
    info = {}
    info["hostname"] = socket.gethostname()
    info["os"] = platform.platform()
    info["python"] = platform.python_version()
    return info
