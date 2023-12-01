#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:20:13
@File: swanlab\env.py
@IDE: vscode
@Description:
    环境变量配置，用于配置一些全局变量，如存储路径等内容
"""
import os
import mimetypes

# 默认存放数据的目录为用户执行python命令时的目录
SWANLAB_FOLDER = os.path.join(os.getcwd(), ".swanlab")
SWANLAB_LOGS_FOLDER = os.path.join(SWANLAB_FOLDER, "logs")
"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
因此，需要指定文件路径与文件名，用于后端服务的响应（这在下面的路由配置中也有说明）
"""
# 注册静态文件路径
mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")
# 静态文件路径
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")
# TODO 后续可以考虑将logo.ico放在assets中，这样就不需要单独响应了
LOGO = os.path.join(TEMPLATE_PATH, "logo.ico")
