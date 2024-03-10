#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-19 15:40:58
@File: swanlab/db/migrate/tag_sort.py
@IDE: vscode
@Description:
    为数据库的Tag表添加sort字段
"""
from peewee import SqliteDatabase, IntegerField
from playhouse.migrate import migrate, SqliteMigrator


def add_sort(db: SqliteDatabase):
    """
    为数据库的Tag表添加sort字段
    """
    migrator = SqliteMigrator(db)
    sort = IntegerField(default=0)
    migrate(migrator.add_column("tag", "sort", sort))
