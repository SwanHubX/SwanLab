#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 17:33:56
@File: swanlab/log/__init__.py
@IDE: vscode
@Description:
    日志记录模块，在设计上swanlog作为一个独立的模块被使用
"""
from .log import SwanLog

swanlog: SwanLog = SwanLog("swanlab")

start_proxy = swanlog.start_proxy

reset = swanlog.reset

__all__ = ["start_proxy", "reset", "swanlog"]
