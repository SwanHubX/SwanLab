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
)


def connect():
    """
    连接数据库，只有调用此方法以后，数据库才会被创建，所有导出的类才可用
    这样设计的原因是因为路径问题，这里需要动态导入settings
    """
    from .settings import swandb

    # 动态绑定数据库
    swandb.connect()
    tables = [
        Project,
        Experiment,
        Tag,
        Chart,
        Namespace,
        Source,
        Display,
    ]
    swandb.bind(tables)
    swandb.create_tables(tables)
