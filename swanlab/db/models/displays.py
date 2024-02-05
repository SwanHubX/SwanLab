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
from ..error import ExistedError, ForeignChartNotExistedError, ForeignNameNotExistedError
from .charts import Chart
from .namespaces import Namespace


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
    chart_id = ForeignKeyField(Chart, backref="displays", on_delete="CASCADE", on_update="CASCADE", null=False)
    """关联的chart，不可为空"""
    namespace_id = ForeignKeyField(Namespace, backref="displays", on_delete="CASCADE", on_update="CASCADE", null=False)
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
        ForeignChartNotExistedError
            外键对应的id不存在(chart不存在)
        ForeignNameNotExistedError
            外键对应的id不存在(namespace不存在)
        ExistedError
           "chart_id" 和 "namespace_id" 的对应关系已存在
        """

        # 检查外键存在性
        if not Chart.filter(Chart.id == chart_id).exists():
            raise ForeignChartNotExistedError("图表不存在")
        if not Namespace.filter(Namespace.id == namespace_id).exists():
            raise ForeignNameNotExistedError("命名空间不存在")

        # 如果sort为None，则自动添加到最后
        if sort is None:
            # 获取当前namespace下的最大索引,如果没有，则为0
            sort = cls.select(fn.Max(cls.sort)).where(cls.namespace_id == namespace_id).scalar()
            sort = 0 if sort is None else sort + 1
        elif sort >= 0 and cls.filter(cls.sort == sort).exists():
            # 先判断当前 sort 是否存在，如果存在，则将大于等于 sort 的索引加1，这样可以直接将数据插入
            if cls.filter(cls.sort == sort).exists():
                cls.update(sort=cls.sort + 1).where(cls.sort >= sort).execute()

        try:
            return super().create(
                chart_id=chart_id,
                namespace_id=namespace_id,
                sort=sort,
                more=cls.dict2json(more),
            )
        except IntegrityError:
            raise ExistedError("display对应关系已存在")
