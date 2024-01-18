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
        # 写入check约束，name和experiment_id加起来唯一
        indexes = ((("tag_id", "chart_id"), True),)

    id = IntegerField(primary_key=True)
    """此表唯一id"""
    tag_id = ForeignKeyField(Tag, backref="sources", null=False)
    """代表关联的tag，不可为空"""
    chart_id = ForeignKeyField(Chart, backref="sources", null=False)
    """关联的chart，不可为空"""
    sort = IntegerField(default=0)
    """此tag在此chart中的排序，值越小越靠前"""
    error = TextField(null=True)
    """代表这个tag数据转化到chart上失败了，错误信息，是一个json字符串"""
    more = TextField(null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新时间"""

    @classmethod
    def create(
        cls,
        tag_id: int,
        chart_id: int,
        error: dict = None,
        more: dict = None,
    ):
        """添加行数据

        Parameters
        ----------
        tag_id : int
            对应的tag_id
        chart_id : int
            对应的chart_id
        error : str
            对应的错误信息
        more : str
            更多信息，json格式
        """
        # TODO 错误检查
        current_time = create_time()
        return super().create(
            tag_id=tag_id,
            chart_id=chart_id,
            error=cls.dict2json(error),
            more=cls.dict2json(more),
            create_time=current_time,
            update_time=current_time,
        )
