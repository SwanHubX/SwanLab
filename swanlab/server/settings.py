#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:43:21
@File: swanlab/server/settings.py
@IDE: vscode
@Description:
    服务部分基础配置
"""
import os
import mimetypes
from ..env import get_swanlog_dir, get_runtime_project

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
"""
# 注册静态文件路径
mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")
# 静态文件路径
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")
# swanlog文件夹路径
SWANLOG_DIR = get_swanlog_dir()
# project.json文件路径
PROJECT_PATH = get_runtime_project()
