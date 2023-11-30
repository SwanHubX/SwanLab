#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:15:42
@File: swanlab\server\database.py
@IDE: vscode
@Description:
    项目模块，创建项目级别数据库库，接下来针对实验级别的数据在此基础上进行操作
"""
from typing import Union
import os
import sqlite3
from ..env import SWANLAB_FOLDER


class SwanProject(object):
    """SwanDataBase类用于创建并连接数据库，实现数据库的增删改查等操作"""

    def __init__(self, project: str, expriment_name: str):
        """初始化数据库连接, 创建数据库，完成数据库的初始化

        Parameters
        ----------
        project : str
            项目名称，用于区分不同的数据库
        name : str
            实验名称
        """
        # 此时必须保证swanlab文件夹存在
        # 数据库名称，格式为{project}.sqlite3
        db_name = project + ".sqlite3"
        # 创建数据库文件，如果已经存在则不创建
        db_path = os.path.join(SWANLAB_FOLDER, db_name)
        # TODO 在此处连接数据库，作为一个属性，接下来
        self.__con = sqlite3.connect(db_path)
        print("db initialized")

    def __create(self):
        """创建数据库"""
        pass

    def add(self, namespace: str, tag: str, index: int, data: Union[int, float]):
        """添加数据到数据库，保存数据

        Parameters
        ----------
        namespace : str
            命名空间，用于区分不同的数据资源
        tag : str
            数据标签，用于区分同一资源下不同的数据
        index : int
            数据索引，用于区分同一资源下同一数据的不同时间点的数据
        data : Union[int, float]
            定位到的数据，暂时只支持int和float类型
        """
        db_path = os.path.join(BASE_PATH, namespace + ".db")
        con = sqlite3.connect(db_path)
        cur = con.cursor()

        table = cur.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{tag}'")
        if table.fetchone() is None:
            cur.execute(f"CREATE TABLE {tag} (id INTEGER PRIMARY KEY AUTOINCREMENT, indes INTEGER, date FLOAT)")
            print(f"Created {tag} table")
        else:
            pass

        # 插入数据
        cur.execute(f"INSERT INTO {tag} VALUES (NULL, {index}, {data})")

        con.commit()
        con.close()

    def query(self, namespace: str, tag: str, index: int):
        """查询数据

        Parameters
        ----------
        namespace : str
            命名空间，用于区分不同的数据资源
        tag : str
            数据标签，用于区分同一资源下不同的数据
        index : int
            数据索引，用于区分同一资源下同一数据的不同时间点的数据
        """
        db_path = os.path.join(BASE_PATH, namespace + ".db")
        con = sqlite3.connect(db_path)
        cur = con.cursor()

        table = cur.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{tag}'")
        if table.fetchone() is None:
            print(f"Table {tag} does not exist")
        else:
            print(f"Table {tag} exists")

        # 查出表中所有数据
        cur.execute(f"SELECT * FROM {tag}")
        res = cur.fetchall()
        print("==================")
        print(res)

        con.close()
