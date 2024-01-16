#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:54:30
@File: swanlab\db\modules\project.py
@IDE: vscode
@Description:
    项目表，对应于0.1.5之前版本中的 project.json
"""
from peewee import CharField, IntegerField
from ..setting import SwanModel


class Project(SwanModel):
    """项目表"""

    id = IntegerField(primary_key=True)
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255, null=True)
    sum = IntegerField(null=True)
    charts = IntegerField(default=0, choices=[0, 1], null=False)
    more = CharField(null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @property
    def get_info(self):
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "sum": self.sum,
            "charts": self.charts,
            "more": self.more,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }
