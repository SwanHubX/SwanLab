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

swanlog: Optional["SwanLog"] = SwanLog("swanlab")

install = swanlog.install

uninstall = swanlog.uninstall

__all__ = ["install", "uninstall", "swanlog"]
