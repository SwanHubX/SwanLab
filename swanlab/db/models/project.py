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

# 默认的项目id应该是1
DEFAULT_PROJECT_ID = 1


class Project(SwanModel):
    """项目表
    在一个工程中，只有一个项目

    Attributes
    ----------
    experiments: list of Experiment
        由 Experiment 表中外键反链接生成的实验列表
    charts: list of Chart
        由 Chart 表中外键反链接生成的图表列表
    namespaces: list of Namespace
        由 Namespace 表中外键反链接生成的命名空间数据列表
    """

    class Meta:
        database = swandb

    id = IntegerField(primary_key=True)
    """项目id"""
    name = CharField(max_length=100, null=False, unique=True)
    """项目名称，唯一"""
    description = CharField(max_length=255, null=True)
    """项目描述，可为空"""
    sum = IntegerField(null=True)
    """项目下实验数量，包括已删除的实验，这是一个只增不减的值"""
    charts = IntegerField(default=0, choices=[0, 1], null=False)
    """是否已经生成项目级别图表，0 未生成，1 已生成"""
    more = CharField(null=True)
    """更多信息配置，json格式"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新时间"""

    @classmethod
    def is_inited(cls):
        """
        判断项目表是否已经初始化，这是为单项目模式设计的
        事实上数据库字段支持多项目，但是目前只需要单项目模式
        """
        return cls.select().count() >= 1

    @classmethod
    def init(cls, name="", description="", more="") -> "Project":
        """初始化项目表，如果已经有项目存在，则不创建，直接返回第一条数据实例

        Parameters
        ----------
        name : str, optional
            项目名称，默认为空字符串
        description : str, optional
            项目描述，默认为空字符串
        more : str, optional
            更多信息配置，默认为空字符串

        Returns
        -------
        Project:
            项目实例
        """
        # 如果已经初始化，则不创建，直接返回第一条数据实例
        if cls.is_inited():
            return cls.filter(cls.id == DEFAULT_PROJECT_ID)[0]
        # 创建项目
        return cls.create(
            name=name,
            description=description,
            more=more,
            create_time=create_time(),
            update_time=create_time(),
        )

    @classmethod
    def increase_sum(cls, id=DEFAULT_PROJECT_ID) -> int:
        """更新实验统计数量，增加1
        Parameters
        ----------
        id : int
            实验id, 默认为DEFAULT_PROJECT_ID

        Returns
        -------
        int:
            当前实验统计数量
        """
        project: "Project" = cls.filter(cls.id == id)[0]
        project.sum += 1
        project.save()
        return project.sum

    @classmethod
    def get_sum(cls, id=DEFAULT_PROJECT_ID) -> int:
        """获取某个项目的历史实验个数"""

        project: "Project" = cls.filter(cls.id == id)[0]
        return project.sum

    @classmethod
    def update_info(cls, name: str, description: str = "", id=DEFAULT_PROJECT_ID) -> "Project":
        """设置实验名和实验描述

        Parameters
        ----------
        name : str
            实验名，不能为空字符串
        description : str
            实验描述，可为空字符串

        Returns
        -------
        Project:
            项目实例
        """
        if name is None or name == "":
            raise ValueError("Invalid project name")
        project: "Project" = cls.filter(cls.id == id)[0]
        project.name = name
        project.description = description
        # 更新更新时间
        project.update_time = create_time()
        project.save()
        return project

    @classmethod
    @SwanModel.result_to_dict
    def get_experiments(cls):
        """获取项目下的所有实验"""

        return cls.select()[0].experiments
