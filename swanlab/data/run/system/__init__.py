#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-29 19:22:55
@File: swanlab\system\__init__.py
@IDE: vscode
@Description:
    硬件监控模块与系统信息获取
"""

from .monitor import SwanSystemMonitor
from .info import get_system_info, get_requirements
