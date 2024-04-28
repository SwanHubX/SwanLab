#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 17:33:56
@File: swanlab/log/__init__.py
@IDE: vscode
@Description:
    日志记录模块，在设计上swanlog作为一个独立的模块被使用，你可以在除了utils的任何地方使用它
"""
from typing import Optional
from .log import SwanLog

swanlog: Optional["SwanLog"] = SwanLog("swanlab")

install = swanlog.install

uninstall = swanlog.uninstall

__all__ = ["install", "uninstall", "swanlog"]
