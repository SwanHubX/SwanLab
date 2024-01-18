#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 20:36:01
@File: swanlab\db\models\sources.py
@IDE: vscode
@Description:
    中间表，多chart多tag
"""

from ..settings import swandb
from .tags import Tag
from .charts import Chart
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField, Check
from ...utils.time import create_time


class Source(SwanModel):
    """多chart多tag"""

    class Meta:
        database = swandb
        indexes = ((("tag_id", "chart_id"), True),)

    id = IntegerField(primary_key=True)
    tag_id = ForeignKeyField(Tag, backref="sources", on_delete="SET NULL", null=True)
    chart_id = ForeignKeyField(Chart, backref="sources", on_delete="SET NULL", null=True)
    error = TextField(null=True)
    more = TextField(null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create(
        cls,
        tag_id: int = None,
        chart_id: int = None,
        error: str = None,
        more: str = None,
    ):
        """添加行数据

        Parameters
        ----------
        tag_id : int, optional
            _description_, by default None
        chart_id : int, optional
            _description_, by default None
        error : str, optional
            _description_, by default None
        more : str, optional
            _description_, by default None
        """

        current_time = create_time()
        return super().create(
            tag_id=tag_id,
            chart_id=chart_id,
            error=error,
            more=more,
            create_time=current_time,
            update_time=current_time,
        )
