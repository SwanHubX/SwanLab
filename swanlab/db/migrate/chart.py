#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-09 17:39:12
@File: swanlab/db/migrate/chart.py
@IDE: vscode
@Description:
    chart相关迁移操作
"""
from peewee import SqliteDatabase, IntegerField
from playhouse.migrate import migrate, SqliteMigrator
from .experiment import (
    add_hidden_opened as add_hidden_opened_experiment,
    add_pinned_opened as add_pinned_opened_experiment,
)
from .project import add_hidden_opened as add_hidden_opened_project, add_pinned_opened as add_pinned_opened_project


def add_status(db):
    """
    为数据库的Chart表添加status字段，默认为0
    """
    migrator = SqliteMigrator(db)
    status = IntegerField(default=0)
    migrate(migrator.add_column("chart", "status", status))
    # 同时，添加sort字段
    add_sort(db)
    # 同时为experiment和project添加hidden_opened、pinned_opened字段
    add_hidden_opened_experiment(db)
    add_pinned_opened_experiment(db)
    add_hidden_opened_project(db)
    add_pinned_opened_project(db)


def add_sort(db):
    """
    为数据库的Chart表添加sort字段，默认为NULL
    """
    migrator = SqliteMigrator(db)
    sort = IntegerField(default=None, null=True)
    migrate(migrator.add_column("chart", "sort", sort))
