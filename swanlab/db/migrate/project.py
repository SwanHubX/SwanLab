#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-09 17:43:27
@File: swanlab/db/migrate/project.py
@IDE: vscode
@Description:
    为数据库的Project表添加迁移操作
"""
from peewee import IntegerField
from playhouse.migrate import migrate, SqliteMigrator


def add_pinned_opened(db):
    """
    为数据库的Project表添加pinned_opened字段,用于记录项目置顶部分是否打开
    默认值为1，表示打开
    """
    migrator = SqliteMigrator(db)
    pinned_opened = IntegerField(default=1, choices=[0, 1])
    migrate(migrator.add_column("project", "pinned_opened", pinned_opened))


def add_hidden_opened(db):
    """
    为数据库的Project表添加hidden_opened字段,用于记录项目隐藏部分是否打开
    默认值为0,表示关闭
    """
    migrator = SqliteMigrator(db)
    hidden_opened = IntegerField(default=0, choices=[0, 1])
    migrate(migrator.add_column("project", "hidden_opened", hidden_opened))
