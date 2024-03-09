#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-05 19:42:52
@File: swanlab/db/migrate/namespace_opened.py
@IDE: vscode
@Description:
    为namespace增加opened字段
"""
from peewee import SqliteDatabase, IntegerField
from playhouse.migrate import migrate, SqliteMigrator


def add_opened(db):
    """
    为数据库的namespace表添加opened字段,用于记录命名空间是否已经被打开
    默认值为False
    """
    migrator = SqliteMigrator(db)
    opened = IntegerField(default=1, choices=[0, 1])
    migrate(migrator.add_column("namespace", "opened", opened))
