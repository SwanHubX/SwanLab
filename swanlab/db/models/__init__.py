#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 13:08:48
@File: swanlab/db/models/__init__.py
@IDE: vscode
@Description:
    在此处自定义模型类，继承自SwanModel
    注意，在数据表模型类中，自实现的方法的命名需要有一定的区分度
    例如 project 表中，删除方法可以命名为 delete_project
    以避免和 peewee 提供的接口冲突
"""
from .charts import Chart
from .displays import Display
from .experiments import Experiment
from .namespaces import Namespace
from .projects import Project
from .tags import Tag
from .sources import Source
