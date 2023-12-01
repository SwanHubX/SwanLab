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

# 默认存放数据的目录为用户执行python命令时的目录
SWANLAB_FOLDER = os.path.join(os.getcwd(), ".swanlab", "logs")
