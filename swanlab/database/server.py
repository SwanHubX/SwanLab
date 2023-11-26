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


class SwanDataBase(object):
    """SwanDataBase类用于创建并连接数据库，实现数据库的增删改查等操作"""

    def __init__(self):
        """初始化数据库连接"""
        self.tmp = {"default": {}}  # 临时变量，用于存储数据库连接，先用字典代替
        pass

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
        # FIXME: 未完成，需要实现数据库连接，目前简单以字典代替
        # 如果命名空间不存在，报错
        if namespace not in self.tmp:
            raise Exception("namespace {} not exist".format(namespace))
        # 如果标签不存在，创建标签
        if tag not in self.tmp[namespace]:
            self.tmp[namespace][tag] = []
        # index在此处没有作用（暂时）
        self.tmp[namespace][tag]: list.append(data)
