#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:15:42
@File: swanlab\server\database.py
@IDE: vscode
@Description:
    数据库模块，用于创建数据库连接并执行一些数据库操作
"""
from typing import Union
import os
import sqlite3

# 默认存放数据的目录
DB_FOLDER = "db"
BASE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DB_FOLDER)

class SwanDataBase(object):
    """SwanDataBase类用于创建并连接数据库，实现数据库的增删改查等操作"""

    def __init__(self):
        """初始化数据库连接"""
        if not os.path.exists(BASE_PATH):
            try:
                os.mkdir(BASE_PATH)
                print(f"目录 '{BASE_PATH}' 创建成功。")
            except FileExistsError:
                pass
            except OSError as e:
                print(f"无法创建目录 '{BASE_PATH}': {e}")
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
        db_path = os.path.join(BASE_PATH, namespace + '.db')
        con = sqlite3.connect(db_path)
        cur = con.cursor()

        table = cur.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{tag}'")
        if table.fetchone() is None:
            cur.execute(
                f"CREATE TABLE {tag} (id INTEGER PRIMARY KEY AUTOINCREMENT, indes INTEGER, date FLOAT)")
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
        db_path = os.path.join(BASE_PATH, namespace + '.db')
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
