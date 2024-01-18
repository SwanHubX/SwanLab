#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 20:43:32
@File: swanlab\db\models\displays.py
@IDE: vscode
@Description:
    中间表，namespace 和 chart，设计排序和排序可变
"""

from ..settings import swandb
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField, Check
from .charts import Chart
from .namespaces import Namespace
from ...utils.time import create_time


class Display(SwanModel):
    """namespace 和 chart，设计排序和排序可变"""

    class Meta:
        database = swandb

    id = IntegerField(primary_key=True)
    chart_id = ForeignKeyField(Chart, backref="displays", on_delete="SET NULL", null=True)
    namespace_id = ForeignKeyField(Namespace, backref="displays", on_delete="SET NULL", null=True)
    index = IntegerField()
    more = TextField(default=None, null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create(
        cls,
        chart_id: int = None,
        namespace_id: int = None,
        index: int = -1,
        more: str = None,
    ):
        """添加行数据

        Parameters
        ----------
        chart_id : int, optional
            _description_, by default None
        namespace_id : int, optional
            _description_, by default None
        index : int, optional
            _description_, by default -1
        more : str, optional
            _description_, by default None

        Returns
        -------
        _type_
            _description_
        """

        current_time = create_time()
        return super().create(
            chart_id=chart_id,
            namespace_id=namespace_id,
            index=index,
            more=more,
            create_time=current_time,
            update_time=current_time,
        )
