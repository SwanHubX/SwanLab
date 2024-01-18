#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 13:52:35
@File: swanlab\db\models\chart.py
@IDE: vscode
@Description:
    实验图标表
"""

from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField
from .experiment import Experiment
from .project import Project
from ...utils.time import create_time


class Chart(SwanModel):
    """chart表"""

    id = IntegerField(primary_key=True)
    experiment_id = ForeignKeyField(Experiment, backref="charts", on_delete="SET NULL")
    project_id = ForeignKeyField(Project, backref="charts", default=1, on_delete="SET NULL")
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255)
    system = IntegerField(default=1, choices=[-1, 0, 1])
    type = CharField(max_length=10, null=False)
    reference = CharField(max_length=10, default="step")
    config = TextField(default="")
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create_chart(
        cls,
        experiment_id,
        name,
        project_id=1,
        description="",
        system=1,
        type="default",
        reference="step",
        config="",
        more="",
    ):
        """创建实验图标
        必传：
            - epxeriment_id：实验id
            - name：图标名称
        """

        try:
            cls.create(
                experiment_id=experiment_id,
                project_id=project_id,
                name=name,
                description=description,
                system=system,
                type=type,
                reference=reference,
                config=config,
                more=more,
                create_time=create_time(),
                update_time=create_time(),
            )
        except Exception as e:
            raise e
