#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 10:52:10
@File: swanlab\db\__init__.py
@IDE: vscode
@Description:
    数据库模型模块，在设计上应该在data模块或者server模块被初始化完毕后导入
    导入之前必须确保文件路径存在，否则会报错
"""
from .models import (
    Project,
    Experiment,
    Tag,
    Chart,
    Namespace,
    Source,
    Display,
)
from .error import (
    ExistedError,
    NotExistedError,
    ForeignProNotExistedError,
    ForeignExpNotExistedError,
    ForeignTagNotExistedError,
    ForeignChartNotExistedError,
    ForeignNameNotExistedError,
    ChartTypeError,
)
from .table_config import tables
from .db_connect import connect
from .utils import add_multi_chart
