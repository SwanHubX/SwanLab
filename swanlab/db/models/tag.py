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
from peewee import IntegrityError
from ..error import ExistedError, NotExistedError
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
    # backref是反向引用，可以通过实验获取tag，比如experiment.tags获取此实验下的所有tag
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
        NotExistedError
            实验对应的id不存在
        ExistedError
            此tag已经在数据库中存在
        """
        # 如果实验id不存在，则抛出异常
        if not Experiment.get_experiment(experiment_id):
            raise NotExistedError("Experiment id not found: {}".format(experiment_id))

        # 尝试创建实验tag，如果已经存在则抛出异常
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
        except IntegrityError:
            raise ExistedError("Tag already exists: {}".format(name))

    @classmethod
    @SwanModel.result_to_list
    def get_tags(cls, experiment_id: int) -> list:
        """获取指定实验下的所有tag
        Parameters
        ----------
        experiment_id : int
            实验id

        Returns
        -------
        list
            返回tag数据列表,每个元素是一个字典, 代表一条tag数据
            如果对应实验不存在tag，则返回空列表

        Raises
        ------
        NotExistedError
            实验对应的id不存在
        """
        if not Experiment.get_experiment(experiment_id):
            raise NotExistedError("Experiment id not found: {}".format(experiment_id))
        return cls.filter(cls.experiment_id == experiment_id)

    @classmethod
    @SwanModel.result_to_dict
    def get_tag(cls, id):
        """根据id获取tag
        Parameters
        ----------
        id : int
            tag的id
        """
        try:
            return cls.filter(cls.id == id)
        except IndexError:
            raise NotExistedError("Tag id not found: {}".format(id))

    @classmethod
    @SwanModel.result_to_dict
    def get_tag_by_name(cls, experiment_id: int, name: str):
        """
        已知实验id和tag名称，获取tag数据
        Parameters
        ----------
        experiment_id : int
            实验id
        name : str
            tag名称

        """

        return cls.filter(cls.experiment_id == experiment_id, cls.name == name)
