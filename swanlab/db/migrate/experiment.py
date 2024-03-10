#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-29 14:34:29
@File: swanlab/db/migrate/exp_finish_time.py
@IDE: vscode
@Description:
    为数据库的Experiment表添加finish_time字段，用于记录实验结束时间
    实验结束时间是指实验被标记为完成/崩溃的时间
"""

from peewee import SqliteDatabase, CharField, IntegerField
from playhouse.migrate import migrate, SqliteMigrator


def add_finish_time(db):
    """
    为数据库的Experoment表添加finish_time字段,用于记录实验结束时间
    默认值为当前行的update_time
    """
    migrator = SqliteMigrator(db)
    finish_time = CharField(max_length=30, null=True, default=None)
    migrate(migrator.add_column("experiment", "finish_time", finish_time))
    # 原子化操作，遍历所有status不为0的实验，将其finish_time设置为update_time
    for exp in db.execute_sql("select id, update_time from experiment where status != 0"):
        db.execute_sql(f"update experiment set finish_time = '{exp[1]}' where id = {exp[0]}")


def add_pinned_opened(db):
    """
    为数据库的Experiment表添加pinned_opened字段,用于记录实验置顶部分是否打开
    默认值为1，表示打开
    """
    migrator = SqliteMigrator(db)
    pinned_opened = IntegerField(default=1, choices=[0, 1])
    migrate(migrator.add_column("experiment", "pinned_opened", pinned_opened))


def add_hidden_opened(db):
    """
    为数据库的Experiment表添加hidden_opened字段,用于记录实验隐藏部分是否打开
    默认值为0,表示关闭
    """
    migrator = SqliteMigrator(db)
    hidden_opened = IntegerField(default=0, choices=[0, 1])
    migrate(migrator.add_column("experiment", "hidden_opened", hidden_opened))
