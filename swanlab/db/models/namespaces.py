#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 15:12:37
@File: swanlab\db\models\namespace.py
@IDE: vscode
@Description:
    实验表格命名空间：一个实验/项目下可以有多个命名空间
"""
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField, IntegrityError, Check, fn, DatabaseProxy
from ..error import ExistedError, ForeignProNotExistedError, ForeignExpNotExistedError
from .experiments import Experiment
from .projects import Project


class Namespace(SwanModel):
    """命名空间表

    Parameters
    ----------
    SwanModel : Class
        自建基类
    """

    class Meta:
        database = DatabaseProxy()

        """
        name 和 experiment_id 和 project_id 组合后需要确保唯一性
        """

        # 同一个项目/实验下，命名空间名称不能重复
        indexes = ((("name", "experiment_id"), True), (("name", "project_id"), True))

        constraints = [
            # project_id 和 experiment_id 不能同时为空，也不能同时不为空
            Check(
                "(project_id IS NULL AND experiment_id IS NOT NULL)"
                " OR "
                "(project_id IS NOT NULL AND experiment_id IS NULL)"
            ),
            # index必须大于等于0
            Check("sort >= 0"),
        ]

    id = IntegerField(primary_key=True)
    """namespace唯一id"""
    experiment_id = ForeignKeyField(
        Experiment, backref="namespaces", on_delete="CASCADE", on_update="CASCADE", null=True
    )
    """外键，对应的实验id，与project_id只能有一个为null"""
    project_id = ForeignKeyField(Project, backref="namespaces", on_delete="CASCADE", on_update="CASCADE", null=True)
    """外键，对应的项目id，与experiment_id只能有一个为null"""
    name = CharField(max_length=100, null=False)
    """这个命名空间的名称，同一个项目/实验下，命名空间名称不能重复"""
    description = CharField(max_length=255, null=True)
    """命名空间描述，可为空"""
    sort = IntegerField()
    """命名空间索引，用于排序，同一个项目/实验下，命名空间索引不能重复，索引越小，排序越靠前，索引>=0"""
    opened = IntegerField(default=1, choices=[0, 1])
    """命名空间是否已经被打开，默认为1，表示已经被打开，0表示未被打开"""
    more = TextField(default=None, null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新的时间"""

    def __dict__(self):
        return {
            "id": self.id,
            "experiment_id": self.experiment_id.__dict__(),
            "project_id": self.project_id,
            "name": self.name,
            "description": self.description,
            "sort": self.sort,
            "more": self.more,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

    @classmethod
    def create(
        cls,
        name: str,
        experiment_id: int = None,
        project_id: int = None,
        description: str = None,
        sort=None,
        more: dict = None,
    ) -> "Namespace":
        """创建命名空间

        Parameters
        ----------
        name : str
            命名空间名称
        experiment_id : int, optional
            实验id, 默认为None，但是实验id和项目id不能同时为None，也不能同时不为None
        project_id : int, optional
            项目id, 默认为None，但是实验id和项目id不能同时为None，也不能同时不为None
        description : str, optional
            命名空间描述, 默认为""
        sort : int, optional
            命名空间索引, >=0, 越小排在越前面，默认为None，自动添加到最后
            如果指定了sort，那么将会把其他sort大于等于指定值的命名空间索引+1

        Returns
        -------
        Namespace : Namespace
            创建的命名空间

        Raises
        ------
        ForeignProNotExistedError
            外键对应的id不存在(project不存在)
        ForeignExpNotExistedError
            外键对应的id不存在(experiment不存在)
        ExistedError
           同名命名空间不可存在于同一个实验/项目中
        """

        # 检查外键存在性
        if experiment_id and not Experiment.filter(Experiment.id == experiment_id).exists():
            raise ForeignExpNotExistedError("实验不存在")
        if project_id and not Project.filter(Project.id == project_id).exists():
            raise ForeignProNotExistedError("项目不存在")

        if sort is None:
            # 如果没有指定sort，那么就自动添加到最后
            # 寻找当前实验/项目下的最大索引，如果实验存在就在实验下找，如果实验不存在就在项目下找
            if experiment_id is not None:
                sort = cls.select(fn.Max(cls.sort)).where(cls.experiment_id == experiment_id).scalar()
            else:
                sort = cls.select(fn.Max(cls.sort)).where(cls.project_id == project_id).scalar()
            sort = 0 if sort is None else sort + 1
        elif sort >= 0:
            # 如果指定了sort，那么将会把其他sort大于等于指定值的命名空间索引+1
            if experiment_id is not None:

                cls.update(sort=cls.sort + 1).where(
                    cls.experiment_id == experiment_id,
                    cls.name != name,
                    cls.sort >= sort,
                ).execute()
            else:
                cls.update(sort=cls.sort + 1).where(
                    cls.project_id == project_id,
                    cls.name != name,
                    cls.sort >= sort,
                ).execute()

        try:
            return super().create(
                name=name,
                experiment_id=experiment_id,
                project_id=project_id,
                description=description,
                sort=sort,
                more=more,
            )
        except IntegrityError:
            raise ExistedError("命名空间已存在")
