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
from peewee import SqliteDatabase
from .error import (
    ExistedError,
    NotExistedError,
)
from ..env import get_db_path
import os


tables = [Project, Experiment, Tag, Chart, Namespace, Source, Display]

# 判断是否已经binded了
binded = False


def connect(autocreate=False) -> SqliteDatabase:
    """
    连接数据库，只有调用此方法以后，数据库才会被创建，所有导出的类才可用
    这样设计的原因是因为路径问题，这里需要动态导入settings
    第一次调用此方法时会创建数据库，返回None，以后调用此方法会返回数据库实例，通过其自身的`connect`方法连接数据库
    完成复杂的数据库操作后，需要调用其`close`方法关闭数据库连接

    Parameters:
    ----------
    existed : bool
        是否指定创建数据库，如果为True，则会自动创建数据库，否则不会自动创建数据库
        设置为False时，如果目标数据库不存在，则会抛出FileNotFoundError异常
        但是即使设置为true，文件夹不存在的情况下也会抛出FileNotFoundError异常

    Return:
    -------
    swandb :
        数据库实例

    Raises:
    -------
    FileNotFoundError :
        如果数据库不存在且existed为False，则会抛出FileNotFoundError异常
        或者路径不存在，也会抛出FileNotFoundError异常
    """
    global binded
    path = get_db_path()
    # 判断数据库是否存在
    db_exists = os.path.exists(path)
    path_exists = os.path.exists(os.path.dirname(path))
    if not path_exists or (not db_exists and not autocreate):
        raise FileNotFoundError(f"DB file {path} not found")

    swandb = SqliteDatabase(path, pragmas={"foreign_keys": 1})
    if not binded:
        # 动态绑定数据库
        swandb.connect()
        swandb.bind(tables)
        swandb.create_tables(tables)
        swandb.close()
        binded = True
    return swandb
