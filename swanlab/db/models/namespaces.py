#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 15:12:37
@File: swanlab\db\models\namespace.py
@IDE: vscode
@Description:
    实验表格命名空间：一个实验/项目下可以有多个命名空间
"""

from ..settings import swandb
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField
from .experiments import Experiment
from .projects import Project


class Namespace(SwanModel):
    """命名空间表

    Parameters
    ----------
    SwanModel : Class
        自建基类
    """

    class Meta:
        database = swandb
        # name 和 experiment_id 和 project_id 组合后唯一
        indexes = (("name", "experiment_id", "project_id"), True)

    id = IntegerField(primary_key=True)
    experiment_id = ForeignKeyField(Experiment, backref="namespaces", on_delete="SET NULL")
    project_id = ForeignKeyField(Project, backref="namespaces", on_delete="SET NULL")
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255, null=True)
    index = IntegerField()
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)
