#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-18 21:28:47
@File: swanlab/data/auth/__init__.py
@IDE: vscode
@Description:
    用户认证模块，当用户登录时，需要验证用户的身份
    这一块因为可能swanlog模块没有初始化，所以需要自己单独打印一下
"""
from .login import terminal_login, code_login
