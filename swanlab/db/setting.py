#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 11:00:55
@File: swanlab\db\setting.py
@IDE: vscode
@Description:
    一些配置
"""
from peewee import SqliteDatabase, Model

database_name = "swanlab.db"

# 连接到SQLite数据库，设置数据库文件名为swanlab.db
swandb = SqliteDatabase(database_name)


class SwanModel(Model):
    """基础模型类，用于定义数据库表的基本信息"""

    class Meta:
        """设置元信息"""

        database = swandb
