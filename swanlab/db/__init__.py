#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:52:10
@File: swanlab\db\__init__.py
@IDE: vscode
@Description:
    数据库模型模块，在设计上应该在data模块或者server模块被初始化完毕后导入
"""
from .setting import swandb
from .modules.project import Project
from .modules.experiment import Experiment
from .modules.tag import Tag


swandb.connect()
# 在连接时提前创建，确保在查询时数据表一定存在
swandb.create_tables(
    [
        Project,
        Experiment,
        Tag,
    ],
    safe=True,
)
