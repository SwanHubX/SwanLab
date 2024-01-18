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
from peewee import CharField, IntegerField, ForeignKeyField, TextField
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
        # name 和 experiment_id 和 project_id 组合后唯一
        # indexes = ((("name", "experiment_id", "project_id"), True),)

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
    def create_namespace(
        cls,
        name,
        experiment_id=None,
        project_id=None,
        description="",
        more="",
        index=1,
    ):
        """创建命名空间

        Parameters
        ----------
        name : str
            命名空间名称
        experiment_id : int, optional
            实验ID
        project_id : int, optional
            项目ID
        description : str, optional
            命名空间描述
        more : str, optional
            额外信息
        index : int, optional
            命名空间索引

        Returns
        -------
        Namespace : Namespace
            创建的命名空间
        """

        current_time = create_time()
        kwargs = {
            "experiment_id": 1,
            "project_id": None,
            "name": name,
            "description": description,
            "more": more,
            "index": index,
            "create_time": current_time,
            "update_time": current_time,
        }

        # if project_id is not None:
        #     kwargs["project_id"] = project_id
        # elif experiment_id is not None:
        #     kwargs["experiment_id"] = experiment_id
        # else:
        #     raise ValueError("experiment_id and project_id is None")

        return cls.create(
            name="name",
            index=1,
            create_time=current_time,
            update_time=current_time,
            experiment_id=None,
            project_id=None,
        )
