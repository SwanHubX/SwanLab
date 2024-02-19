#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 20:43:22
@File: swanlab\db\modules\tag.py
@IDE: vscode
@Description:
    实验 tag 表
"""
from ..model import SwanModel
from peewee import ForeignKeyField, CharField, TextField, IntegerField, IntegrityError, DatabaseProxy, fn
from ..error import ExistedError, ForeignExpNotExistedError
from .experiments import Experiment


class Tag(SwanModel):
    """实验tag数据表

    Parameters
    ----------
    SwanModel : class
        自封装类
    """

    class Meta:
        database = DatabaseProxy()
        # 写入check约束，name和experiment_id加起来唯一
        indexes = ((("name", "experiment_id"), True),)

    id = IntegerField(primary_key=True)
    """tag唯一id"""
    # backref是反向引用，可以通过实验获取tag，比如experiment.tags获取此实验下的所有tag
    experiment_id = ForeignKeyField(Experiment, backref="tags", null=False, on_delete="CASCADE", on_update="CASCADE")
    """tag对应的实验id，代表这个tag由谁创建"""
    name = CharField(max_length=255, null=False)
    """tag名称，同一个实验下，tag名称不能重复"""
    type = CharField(max_length=10, null=False)
    """tag的类型"""
    description = CharField(max_length=100, null=True)
    """tag的描述，可为空"""
    system = IntegerField(default=0, choices=[0, 1])
    """标识这个tag数据由系统生成还是用户生成，0: 用户生成，1: 系统生成，默认为0"""
    sort = IntegerField(default=0)
    """tag在实验中的排序，值越小越靠前"""
    more = TextField(null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""
    create_time = CharField(max_length=30, null=False)
    """tag的创建时间"""
    update_time = CharField(max_length=30, null=False)
    """tag的更新时间"""

    def __dict__(self):
        return {
            "id": self.id,
            "experiment_id": self.experiment_id,
            "name": self.name,
            "description": self.description,
            "system": self.system,
            "more": self.more,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

    @classmethod
    def create(
        cls,
        experiment_id: int,
        name: str,
        type: str,
        description="",
        system: int = 0,
        more: dict = None,
    ) -> "Tag":
        """在tags表中创建一条tag指定实验的tag数据

        Parameters
        ----------
        experiment_id : int
            这将成为tag数据的外键，指向实验表中的某一条实验数据
        name : str
            tag的名称，字符串类型，小于255个字符
        type : str
            tag的类型，字符串类型，小于10个字符
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
        ForeignExpNotExistedError
            实验对应的id不存在
        ExistedError
            此tag已经在数据库中存在（tag-experiment_id唯一）
        """
        # 如果实验id不存在，则抛出异常
        if not Experiment.filter(Experiment.id == experiment_id).exists():
            raise ForeignExpNotExistedError("experiment不存在")
        # 获取当前实验下tag的最大排序索引，如果没有则为0
        sort = Tag.select(fn.Max(Tag.sort)).where(Tag.experiment_id == experiment_id).scalar()
        sort = sort + 1 if sort is not None else 0
        # 尝试创建实验tag，如果已经存在则抛出异常
        try:
            return super().create(
                experiment_id=experiment_id,
                name=name,
                type=type,
                description=description,
                system=system,
                sort=sort,
                more=more,
            )
        except IntegrityError:
            raise ExistedError("tag已存在")
