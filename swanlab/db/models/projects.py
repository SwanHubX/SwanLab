#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:54:30
@File: swanlab\db\modules\project.py
@IDE: vscode
@Description:
    项目表，对应于0.1.5之前版本中的 project.json
"""
from peewee import CharField, IntegerField, DatabaseProxy
from ..model import SwanModel
from swanlab.package import get_package_version


class Project(SwanModel):
    """项目表
    目前，在一个工程中，只有一个项目
    """

    # 默认的项目id应该是1
    DEFAULT_PROJECT_ID = 1

    class Meta:
        database = DatabaseProxy()

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
    pinned_opened = IntegerField(default=1, choices=[0, 1])
    """多实验图表置顶部分是否打开，默认值为1，表示打开"""
    hidden_opened = IntegerField(default=0, choices=[0, 1])
    """多实验图表隐藏部分是否打开，默认值为0,表示关闭"""
    more = CharField(null=True)
    """更多信息配置，json格式"""
    version = CharField(max_length=30, null=False)
    """创建时的版本号"""
    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新时间"""

    def __dict__(self):
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "sum": self.sum,
            "charts": self.charts,
            "more": self.more,
            "pinned_opened": self.pinned_opened,
            "hidden_opened": self.hidden_opened,
            "version": self.version,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

    @classmethod
    def init(cls, name="", description=None, more="") -> "Project":
        """
        静态方法
        初始化项目表，如果已经有项目存在，则不创建，直接返回第一条数据实例

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
        if cls.select().count() >= 1:
            return cls.filter(cls.id == cls.DEFAULT_PROJECT_ID)[0]
        # 创建项目
        return cls.create(
            name=name,
            description=description,
            more=more,
            version=get_package_version(),
        )

    @classmethod
    def increase_sum(cls, id=DEFAULT_PROJECT_ID) -> int:
        """
        静态方法
        更新实验统计数量，增加1
        此方法通常在创建实验时被Experiment类调用
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
        project.sum = project.sum + 1 if project.sum else 1
        project.save()
        return project.sum

    @classmethod
    def decrease_sum(cls, id=DEFAULT_PROJECT_ID) -> int:
        """
        静态方法
        更新实验统计数量，减少1
        此方法通常在删除实验时被Experiment类调用
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
        project.sum = project.sum - 1 if project.sum else 0
        project.save()
        return project.sum
