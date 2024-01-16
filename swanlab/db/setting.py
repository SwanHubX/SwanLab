#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 11:00:55
@File: swanlab\db\setting.py
@IDE: vscode
@Description:
    一些配置
"""
from peewee import SqliteDatabase

# 连接到SQLite数据库，设置数据库文件名为swanlab.db
swandb = SqliteDatabase("swanlab.db")
