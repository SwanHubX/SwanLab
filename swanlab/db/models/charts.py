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
from peewee import CharField, IntegerField, ForeignKeyField, TextField, IntegrityError, Check, DatabaseProxy, fn
from ..error import ExistedError, ForeignExpNotExistedError, ForeignProNotExistedError
from .experiments import Experiment
from .projects import Project


class Chart(SwanModel):
    """chart表，图表在namspace下的排序交由display表处理"""

    class Meta:
        database = DatabaseProxy()
        # 通过meta规定name和project_id的唯一性
        indexes = ((("name", "experiment_id"), True), (("name", "project_id"), True))

        constraints = [
            # project_id 和 experiment_id 不能同时为空，也不能同时不为空
            Check(
                "(project_id IS NULL AND experiment_id IS NOT NULL) OR (project_id IS NOT NULL AND experiment_id IS NULL)"
            ),
        ]

    id = IntegerField(primary_key=True)
    """图表id, 自增"""
    experiment_id = ForeignKeyField(Experiment, backref="charts", on_delete="CASCADE", on_update="CASCADE", null=True)
    """外键，关联的实验id，与project_id只有一个为NULL"""
    project_id = ForeignKeyField(
        Project, backref="charts", default=1, on_delete="CASCADE", on_update="CASCADE", null=True
    )
    """外键，关联的项目id，与experiment_id只有一个为NULL"""
    name = CharField(max_length=100, null=False)
    """图表名称"""
    description = CharField(max_length=255, null=True)
    """图表描述，可为空"""
    system = IntegerField(default=1, choices=[-1, 0, 1])
    """是否为创建tag时自动生成的图表，-1: 删除的自动生成的图表，0: 否，1: 是，系统图表不可删除，只能改为-1"""
    type = CharField(max_length=10, null=False)
    """图表类型，由创建者决定，这与数据库本身无关"""
    reference = CharField(max_length=10, default="step")
    """图表数据参考，由创建者决定，这与数据库本身无关"""
    status = IntegerField(null=False, default=0)
    """图表状态 0代表默认状态 1代表被pinned -1代表被hidden"""
    sort = IntegerField(null=True, default=None)
    """图表在被pinned或hidden时的排序，越小越靠前，越大越靠后，如果为NULL则代表未被pinned或hidden"""
    config = TextField(null=True)
    """图表的其他配置，这实际上是一个json字符串"""
    more = TextField(null=True)
    """更多信息配置，json格式"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新时间"""

    def __dict__(self):
        experiment_id = None if self.experiment_id is None else self.experiment_id.__dict__()
        project_id = None if self.project_id is None else self.project_id.__dict__()
        return {
            "id": self.id,
            "experiment_id": experiment_id,
            "project_id": project_id,
            "name": self.name,
            "description": self.description,
            "system": self.system,
            "type": self.type,
            "reference": self.reference,
            "config": self.config,
            "status": self.status,
            "sort": self.sort,
            "more": self.more,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

    @classmethod
    def create(
        cls,
        name: str,
        description: str = None,
        experiment_id: int = None,
        project_id: int = None,
        system: int = 1,
        type: str = "default",
        reference: str = "step",
        config: dict = None,
        more: dict = None,
    ) -> "Chart":
        """创建实验图表

        Parameters
        ----------
        name : str
            实验图表名称
        description : str
            实验图表描述
        experiment_id : int, optional
            外键，关联的实验id，与project_id只有一个为NULL, 默认为None，但与project_id只有一个为None
        project_id : int, optional
            外键，关联的项目id，与experiment_id只有一个为NULL, 默认为None，但与experiment_id只有一个为None
        system : int, optional
            是否为创建tag时自动生成的图表，-1: 删除的自动生成的图表，0: 否，1: 是，系统图表不可删除，只能改为-1, 默认为是
        type : str, optional
            图表类型，由创建者决定，这与数据库本身无关, 默认为"default"
        reference : str, optional
            图表数据参考，由创建者决定，这与数据库本身无关, 默认为"step"
        config : dict, optional
            图表的其他配置，这实际上是一个json字符串, 默认为None
        more : dict, optional
            更多信息配置，json格式, 默认为None

        Returns
        -------
        Chart : Chart
            创建的实验图表

        Raises
        -------
        ForeignExpNotExistedError
            实验不存在
        ForeignProNotExistedError
            项目不存在
        ExistedError
            图标必须唯一
        """

        # 检查外键存在性
        if project_id and not Project.filter(Project.id == project_id).exists():
            raise ForeignProNotExistedError("项目不存在")
        if experiment_id and not Experiment.filter(Experiment.id == experiment_id).exists():
            raise ForeignExpNotExistedError("实验不存在")

        try:
            return super().create(
                experiment_id=experiment_id,
                project_id=project_id,
                name=name,
                description=description,
                system=system,
                type=type,
                reference=reference,
                config=config,
                more=more,
            )
        except IntegrityError:
            raise ExistedError("图表已存在")

    @classmethod
    def pin(cls, id: int, sort: int = None) -> "Chart":
        """
        把某个chart置顶

        Parameters
        ----------
        id : int
            图表id
        sort : int
            图表在被pinned或hidden时的排序，越小越靠前，越大越靠后，如果为NULL则代表未被pinned或hidden, 默认为None
            默认排在最后，也就是传入None的时候

        Returns
        -------
        Chart : Chart
            置顶的图表
        """
        chart: Chart = cls.get(cls.id == id)
        chart.status = 1
        if sort is None:
            # 判断当前实验/项目下有多少个hidden的chart
            project_id = None if chart.project_id is None else chart.project_id.id
            experiment_id = None if chart.experiment_id is None else chart.experiment_id.id
            # 获取当前namespace下的最大索引,如果没有，则为0
            sort = (
                cls.select(fn.Max(cls.sort))
                .where(cls.project_id == project_id, cls.experiment_id == experiment_id, cls.status == 1)
                .scalar()
            )
            sort = 0 if sort is None else sort + 1
        chart.sort = sort
        chart.save()
        return chart

    @classmethod
    def restore(cls, id: int) -> "Chart":
        """
        把某个chart恢复正常状态

        Parameters
        ----------
        id : int
            图表id

        Returns
        -------
        Chart : Chart
            恢复的图表
        """
        chart: Chart = cls.get(cls.id == id)
        chart.status = 0
        chart.sort = None
        chart.save()
        return chart

    @classmethod
    def hide(cls, id: int, sort: int = None) -> "Chart":
        """
        把某个chart隐藏

        Parameters
        ----------
        id : int
            图表id
        sort : int
            图表在被pinned或hidden时的排序，越小越靠前，越大越靠后，如果为NULL则代表未被pinned或hidden, 默认为None
            默认排在最后，也就是传入None的时候

        Returns
        -------
        Chart : Chart
            隐藏的图表
        """
        chart: Chart = cls.get(cls.id == id)
        chart.status = -1
        if sort is None:
            # 判断当前实验/项目下有多少个hidden的chart
            project_id = None if chart.project_id is None else chart.project_id.id
            experiment_id = None if chart.experiment_id is None else chart.experiment_id.id
            # 获取当前namespace下的最大索引,如果没有，则为0
            sort = (
                cls.select(fn.Max(cls.sort))
                .where(cls.project_id == project_id, cls.experiment_id == experiment_id, cls.status == -1)
                .scalar()
            )
            sort = 0 if sort is None else sort + 1
        chart.sort = sort
        chart.save()
        return chart
