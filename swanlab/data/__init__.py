#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 16:08:15
@File: swanlab/data/__init__.py
@IDE: vscode
@Description:
    在此处完成回调注册、swanlog注册，并为外界提供api，提供运行时生成的配置
"""
from .modules import (
    BaseType,
)
from .sdk import (
    init,
    log,
    finish,
    config,
)
