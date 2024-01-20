#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 20:43:32
@File: swanlab\db\models\displays.py
@IDE: vscode
@Description:
    中间表，namespace 和 chart，设计排序和排序可变
"""
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField, IntegrityError, Check, fn, DatabaseProxy
from ..error import ExistedError, NotExistedError
from .charts import Chart
from .namespaces import Namespace
from ...utils.time import create_time


class Display(SwanModel):
    """namespace 和 chart，设计排序和排序可变"""

    class Meta:
        database = DatabaseProxy()
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

    def __dict__(self):
        return {
            "id": self.id,
            "chart_id": self.chart_id,
            "namespace_id": self.namespace_id,
            "sort": self.sort,
            "more": self.more,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

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

        Raises
        ------
        NotExistedError
            外键对应的id不存在(chart/namespace不存在)
        ExistedError
           "chart_id" 和 "namespace_id" 的对应关系已存在
        """

        # 检查外键存在性
        if not Chart.filter(Chart.id == chart_id).exists():
            raise NotExistedError("图表不存在")
        if not Namespace.filter(Namespace.id == namespace_id).exists():
            raise NotExistedError("命名空间不存在")

        # 如果sort为None，则自动添加到最后
        if sort is None:
            max_sort = cls.select(fn.MAX(cls.sort)).scalar()
            if max_sort or max_sort == 0:
                sort = max_sort + 1
            else:
                sort = 0
        elif sort >= 0:
            # 先判断当前 sort 是否存在
            if cls.filter(cls.sort == sort).exists():
                cls.update({cls.sort: cls.sort + 1}).where(cls.sort >= sort).execute()
        else:
            raise ValueError("命名空间索引必须大于等于0")

        current_time = create_time()

        try:
            return super().create(
                chart_id=chart_id,
                namespace_id=namespace_id,
                sort=sort,
                more=cls.dict2json(more),
                create_time=current_time,
                update_time=current_time,
            )
        except IntegrityError:
            raise ExistedError("display对应关系已存在")
