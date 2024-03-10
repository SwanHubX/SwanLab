#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-19 15:41:25
@File: swanlab/db/migrate/__init__.py
@IDE: vscode
@Description:
    数据库迁移模块
"""
from .tag import add_sort
from .experiment import add_finish_time
from .namespace import add_opened
from .chart import add_status
