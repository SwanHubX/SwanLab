#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 20:36:33
@File: swanlab/app/__init__.py
@IDE: vscode
@Description:
        在此处引出swanlab的web服务器框架，名为SwanWeb以及一些神奇的路由配置，以完成在库最外层的函数式调用
"""
from swanlab.env import swc
from swanlab.log import swanlog as swl

# 先初始化配置文件和日志对象
swc.init(swc.getcwd(), "server")
swl.init(swc.output, level="debug")

# 导出app对象
from .router import app
