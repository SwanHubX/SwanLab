#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 20:43:22
@File: swanlab\db\modules\tag.py
@IDE: vscode
@Description:
    实验 tag 表
"""

from ..setting import SwanModel
from peewee import ForeignKeyField, CharField, TextField, IntegerField
from ...utils.time import create_time
from .experiment import Experiment


class Tag(SwanModel):
    """实验tag数据表

    Parameters
    ----------
    SwanModel : class
        自封装类
    """

    id = IntegerField(primary_key=True)
    experiment_id = ForeignKeyField(Experiment, backref="tags", null=False)
    name = CharField(unique=True, max_length=100, null=False)
    description = CharField(max_length=100)
    system = IntegerField(default=0, choices=[0, 1])
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create_tag(
        cls,
        experiment_id,
        name,
        description="",
        system=0,
        more="",
    ):
        """创建实验tag

        Raises
        ------
        ValueError
            所属实验不存在
        """

        # 检查实验是否存在
        if not Experiment.get_experiment(experiment_id):
            raise ValueError("Experiment does not exist")

        try:
            return cls.create(
                experiment_id=experiment_id,
                name=name,
                description=description,
                system=system,
                more=more,
                create_time=create_time(),
                update_time=create_time(),
            )
        except Exception as e:
            raise e

    @classmethod
    @SwanModel.result_to_list
    def get_tags(cls, experiment_id):
        """获取指定实验下的所有tag"""

        return cls.filter(cls.experiment_id == experiment_id)
