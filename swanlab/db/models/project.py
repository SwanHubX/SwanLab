#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:54:30
@File: swanlab\db\modules\project.py
@IDE: vscode
@Description:
    项目表，对应于0.1.5之前版本中的 project.json
"""
from ..settings import swandb
from peewee import CharField, IntegerField
from ..model import SwanModel
from ...utils.time import create_time


class Project(SwanModel):
    """项目表
    在一个工程中，只有一个项目

    Attributes
    ----------
    experiments: list of Experiment
        由 Experiment 表中外键反链接生成的实验列表
    charts: list of Chart
        由 Chart 表中外键反链接生成的图表列表
    """

    class Meta:
        database = swandb

    id = IntegerField(primary_key=True)
    name = CharField(max_length=100, null=False)
    description = CharField(max_length=255, null=True)
    sum = IntegerField(null=True)
    charts = IntegerField(default=0, choices=[0, 1], null=False)
    more = CharField(null=True)
    create_time = CharField(max_length=30, null=False)
    update_time = CharField(max_length=30, null=False)

    @classmethod
    def check_project(cls):
        """检查项目是否存在 => 只存在一个项目
        如果已经有项目存在，返回 true
        """

        return cls.select().count() >= 1

    @classmethod
    def init(cls, name="Sample Project", description="This is a sample project.", sum=0, charts=0, more=""):
        """初始化项目表，如果已经有项目存在，则不创建
        若不满足初始化条件，返回 None
        若初始化成功，返回项目实例，可以在项目实例上获取项目信息
        TODO: 项目名等初始值的设置

        Parameters
        ----------
        name : str, optional
            by default "Sample Project"
        description : str, optional
            by default "This is a sample project."
        sum : int, optional
            by default 0
        charts : int, optional
            by default 0
        more : str, optional
            by default ""

        Returns
        -------
        not create: None
        create: Project instance
        """

        if cls.check_project():
            return None
        # 创建项目
        return cls.create(
            name=name,
            description=description,
            sum=sum,
            charts=charts,
            more=more,
            create_time=create_time(),
            update_time=create_time(),
        )

    @classmethod
    def delete_project(cls):
        """清空项目表"""

        return cls.delete().execute()

    @classmethod
    def update_sum(cls, type="increase"):
        """更新实验统计数量
        TODO: 等实验表建立后，可以从实验表中获取实验数量

        Returns
        -------
        int:
            被操作的行数
        """

        project = cls.select()[0]
        if type == "increase":
            project.sum += 1
        else:
            if project.sum == 0:
                raise ValueError("Experiments number is 0")
            project.sum -= 1
        return project.save()

    @classmethod
    def get_sum(cls):
        """获取实验统计个数"""

        project = cls.select()[0]
        return project.sum

    @classmethod
    def update_info(cls, name: str, description: str = ""):
        """设置实验名、实验描述

        Parameters
        ----------
        name : str
            实验名，不为空
        description : str
            实验描述，可为空

        Returns
        -------
        int:
            被操作的行数
        """
        if name is None or name == "":
            raise ValueError("Invalid project name")
        project = cls.select()[0]
        project.name = name
        project.description = description
        return project.save()

    @classmethod
    def update_updatetime(cls):
        """更新项目更新时间

        Returns
        -------
        int:
            被操作的行数
        """

        return cls.update(update_time=create_time()).where(cls.id == 1).execute()

    @classmethod
    @SwanModel.result_to_dict
    def get_experiments(cls):
        """获取项目下的所有实验"""

        return cls.select()[0].experiments
