#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 14:56:38
@File: swanlab\db\modules\experiment.py
@IDE: vscode
@Description:
    实验数据表
"""

from ..setting import SwanModel
from peewee import Model, SqliteDatabase, ForeignKeyField, CharField, IntegerField, TextField
from .project import Project
from ...utils.time import create_time


# 连接到数据库
db = SqliteDatabase("your_database.db")


# 定义模型类
class Experiment(SwanModel):
    """实验表"""

    id = IntegerField(primary_key=True, unique=True)
    project_id = ForeignKeyField(Project, backref="experiments", default=1)  # 假设存在 Project 模型
    run_id = IntegerField()
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255)
    index = IntegerField()
    status = IntegerField(choices=[-1, 0, 1])
    show = IntegerField(default=1, choices=[0, 1])
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False, default=create_time())
    update_time = CharField(max_length=30, null=False, default=create_time())
