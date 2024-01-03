#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 20:36:33
@File: swanlab/app/__init__.py
@IDE: vscode
@Description:
        在此处引出swanlab的web服务器框架，名为SwanWeb以及一些神奇的路由配置，以完成在库最外层的函数式调用
"""
import os
from ..env import get_swanlog_dir
from ..log import register

register(os.path.join(get_swanlog_dir(), "output.log"))

# 导出app对象
from .router import app
