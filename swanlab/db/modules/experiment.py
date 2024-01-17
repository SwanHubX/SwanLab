#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 14:56:38
@File: swanlab\db\modules\experiment.py
@IDE: vscode
@Description:
    实验数据表
"""

from ..setting import SwanModel
from peewee import ForeignKeyField, CharField, IntegerField, TextField
from .project import Project
from ...utils.time import create_time


# 定义模型类
class Experiment(SwanModel):
    """实验表"""

    id = IntegerField(primary_key=True, unique=True)
    project_id = ForeignKeyField(Project, backref="experiments", default=1)
    run_id = IntegerField(unique=True)
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255)
    index = IntegerField()
    status = IntegerField(choices=[-1, 0, 1])
    show = IntegerField(default=1, choices=[0, 1])
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create_experiment(
        cls,
        run_id,
        project_id=1,
        name="Experiment",
        description="This is a sample experiment.",
        index=0,
        status=0,
        show=1,
        more="",
    ):
        """创建实验表
        其中需要注意：run_id 必须传递！
        TODO: 是否需要在创建实验数据之前 init Project 表？
        """

        # 确保 run_id 唯一
        if not cls.check_run_id(run_id):
            raise ValueError("run_id already exists")
        # 获取之前的实验数量
        experiment_count = cls.select().count()
        experiment = cls.create(
            project_id=project_id,
            run_id=run_id,
            name=name + f"-{experiment_count}",
            description=description,
            index=index,
            status=status,
            show=show,
            more=more,
            create_time=create_time(),
            update_time=create_time(),
        )

        if experiment:
            Project.update_sum(type="increase")
            return experiment
        else:
            raise Exception("Failed to create experiment")

    @classmethod
    def check_run_id(cls, run_id):
        """检查run_id是否可用"""

        return cls.select().where(cls.run_id == run_id).count() == 0

    @classmethod
    @SwanModel.result_to_dict
    def get_experiment(cls, id):
        """根据id获取实验"""

        return cls.select().where(cls.id == id)

    @classmethod
    @SwanModel.result_to_dict
    def get_experiment_by_runid(cls, run_id):
        """根据run_id获取实验"""

        return cls.select().where(cls.run_id == run_id)

    @classmethod
    @SwanModel.result_to_dict
    def get_experiment_by_name(cls, name):
        """根据name获取实验"""

        return cls.select().where(cls.name == name)

    @classmethod
    def delete_experiment(cls, id):
        """删除实验"""

        experiment = cls.select().where(cls.id == id).first()
        if experiment:
            experiment.delete_instance(recursive=True)
            Project.update_sum(type="decrease")
            return True
        else:
            return False

    @classmethod
    def update_info(cls, id: int, name: str, description: str = ""):
        """更新实验信息

        Parameters
        ----------
        id : int
        name : str
        description : str, optional
        """

        experiment = cls.select().where(cls.id == id).first()
        if experiment:
            experiment.name = name
            experiment.description = description
            experiment.save()
            return True
        else:
            raise Exception("Failed to update experiment")

    @classmethod
    def update_updatetime(cls, id: int):
        """实验有更新，刷新 update_time

        Parameters
        ----------
        id : int
        """

        experiment = cls.select().where(cls.id == id).first()
        if experiment:
            experiment.update_time = create_time()
            experiment.save()
            return True
        else:
            raise Exception("Failed to renew update time")

    @classmethod
    def update_status(cls, id: int, status: int):
        """更新实验状态

        Parameters
        ----------
        id : int
        status : int
            * -1: crushed
            *  0: running
            *  1: finished
        """

        experiment = cls.select().where(cls.id == id).first()
        if experiment:
            experiment.status = status
            experiment.save()
            return True
        else:
            raise Exception("Failed to update experiment status")
