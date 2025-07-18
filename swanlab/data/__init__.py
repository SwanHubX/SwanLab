#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 16:08:15
@File: swanlab/data/__init__.py
@IDE: vscode
@Description:
    在此处完成回调注册、swanlog注册，并为外界提供api，提供运行时生成的配置
"""
from .run import (
    SwanLabRun as Run,
    SwanLabRunState as State,
    get_run,
    get_config,
    get_url,
    get_project_url,
    config,
)
from .sdk import (
    login,
    register_callbacks,
    init,
    log,
    finish,
    merge_settings,
)
