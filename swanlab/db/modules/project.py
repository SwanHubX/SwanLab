#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:54:30
@File: swanlab\db\modules\project.py
@IDE: vscode
@Description:
    项目表，对应于0.1.5之前版本中的 project.json
"""
from peewee import Model, CharField, IntegerField
from ..setting import swandb


class Project(Model):
    id = IntegerField(primary_key=True)
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255, null=True)
    sum = IntegerField(null=True)
    charts = IntegerField(default=0, choices=[0, 1], null=False)
    more = CharField(null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    class Meta:
        database = swandb
