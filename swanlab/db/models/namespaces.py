#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 15:12:37
@File: swanlab\db\models\namespace.py
@IDE: vscode
@Description:
    实验表格命名空间：一个实验/项目下可以有多个命名空间
"""

from ..settings import swandb
from ..model import SwanModel
from peewee import CharField, IntegerField, ForeignKeyField, TextField, Check
from .experiments import Experiment
from .projects import Project
from ...utils.time import create_time


class Namespace(SwanModel):
    """命名空间表

    Parameters
    ----------
    SwanModel : Class
        自建基类
    """

    class Meta:
        database = swandb

        """
        name 和 experiment_id 和 project_id 组合后需要确保唯一性
        """
        indexes = (
            # 同一个项目/实验下，命名空间名称不能重复
            (("name", "experiment_id"), True),
            (("name", "project_id"), True),
        )

        constraints = [
            # project_id 和 experiment_id 不能同时为空，也不能同时不为空
            Check(
                "(project_id IS NULL AND experiment_id IS NOT NULL) OR (project_id IS NOT NULL AND experiment_id IS NULL)"
            ),
            # index必须大于等于0
            Check("index >= 0"),
        ]

    id = IntegerField(primary_key=True)
    experiment_id = ForeignKeyField(Experiment, backref="namespaces", on_delete="SET NULL", null=True)
    project_id = ForeignKeyField(Project, backref="namespaces", on_delete="SET NULL", null=True)
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255, null=True)
    index = IntegerField()
    more = TextField(default=None, null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create(
        cls,
        name: str,
        experiment_id: int = None,
        project_id: int = None,
        description: str = None,
        index=None,
    ):
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
        index : int, optional
            命名空间索引, 如果传入的index为-1，则会自动在原本的基础上加1，如果为None，则会自动设置为0

        Returns
        -------
        Namespace : Namespace
            创建的命名空间
        """

        current_time = create_time()
        # TODO 如果index为None，则自动添加到最后

        return super().create(
            name=name,
            experiment_id=experiment_id,
            project_id=project_id,
            description=description,
            index=index,
            create_time=current_time,
            update_time=current_time,
        )
