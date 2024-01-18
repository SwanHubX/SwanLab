#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 20:43:22
@File: swanlab\db\modules\tag.py
@IDE: vscode
@Description:
    实验 tag 表
"""
from ..settings import swandb
from ..model import SwanModel
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

    class Meta:
        database = swandb
        # 写入check约束，name和experiment_id加起来唯一
        indexes = ((("name", "experiment_id"), True),)

    id = IntegerField(primary_key=True)
    experiment_id = ForeignKeyField(Experiment, backref="tags", null=False)
    name = CharField(max_length=255, null=False)
    description = CharField(max_length=100)
    system = IntegerField(default=0, choices=[0, 1])
    more = TextField(default="")
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def create_tag(
        cls,
        experiment_id: int,
        name: str,
        description="",
        system: int = 0,
        more: dict = {},
    ) -> "Tag":
        """在tags表中创建一条tag指定实验的tag数据

        Parameters
        ----------
        experiment_id : int
            这将成为tag数据的外键，指向实验表中的某一条实验数据
        name : str
            tag的名称，字符串类型，小于255个字符
        description : str, optional
            tag的描述，字符串类型，小于100个字符, 默认为空字符串
        system : int, optional
            判断是否是系统自动生成的tag，默认为否
        more : str, optional
            这必须是字典格式数据，代表更多信息, 默认为空字典

        Returns
        -------
        Tag
            返回创建的tag数据实例，可以通过实例获取tag数据

        Raises
        ------
        ValueError
            实验不存在
        ExistedError
            此name对应的
        """

        # 检查实验是否存在
        if not Experiment.get_experiment(experiment_id):
            raise ValueError("Experiment does not exist: {}".format(experiment_id))
        # 尝试创建实验tag
        return cls.create(
            experiment_id=experiment_id,
            name=name,
            description=description,
            system=system,
            more=more,
            create_time=create_time(),
            update_time=create_time(),
        )

    @classmethod
    @SwanModel.result_to_list
    def get_tags(cls, experiment_id):
        """获取指定实验下的所有tag"""

        return cls.filter(cls.experiment_id == experiment_id)

    @classmethod
    @SwanModel.result_to_dict
    def get_tag(cls, id):
        """根据id获取tag"""

        return cls.filter(cls.id == id)

    @classmethod
    @SwanModel.result_to_dict
    def get_tag_by_name(cls, experiment_id, name):
        """根据name获取tag"""

        return cls.filter(cls.experiment_id == experiment_id, cls.name == name)
