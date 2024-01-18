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
        # chart_id和namespace_id加起来唯一
        indexes = ((("chart_id", "namespace_id"), True),)
        # 写入check约束，sort必须大于等于0
        constraints = [Check("sort >= 0")]

    id = IntegerField(primary_key=True)
    """display表唯一id"""
    chart_id = ForeignKeyField(Chart, backref="displays", null=False)
    """关联的chart，不可为空"""
    namespace_id = ForeignKeyField(Namespace, backref="displays", null=False)
    """关联的namespace，不可为空"""
    sort = IntegerField()
    """当前chart在namespace下的排序，索引越小，排序越靠前，索引>=0"""
    more = TextField(default=None, null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新时间"""

    @classmethod
    def create(
        cls,
        chart_id: int,
        namespace_id: int,
        sort: int = None,
        more: dict = None,
    ) -> "Display":
        """添加行数据

        Parameters
        ----------
        chart_id : int
            chart表中的id
        namespace_id : int
            namespace表中的id
        sort : int
            排序索引，索引越小，排序越靠前，索引>=0，如果为None，则自动加到最后
        more : str
            更多信息配置，json格式，将在表函数中检查并解析

        Returns
        -------
        Display
            返回创建的display对象
        """
        current_time = create_time()

        return super().create(
            chart_id=chart_id,
            namespace_id=namespace_id,
            sort=sort,
            more=cls.dict2json(more),
            create_time=current_time,
            update_time=current_time,
        )
