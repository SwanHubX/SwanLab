#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:15:42
@File: swanlab\server\database.py
@IDE: vscode
@Description:
    项目模块，创建项目级别数据库库，接下来针对实验级别的数据在此基础上进行操作
"""
from typing import Union, List
import os
import sqlite3
from ..env import SWANLAB_FOLDER
from .expriments_name import generate_random_tree_name, check_expriment_name, make_expriment_name_unique


class SwanProject(object):
    """SwanDataBase类用于创建并连接数据库，实现数据库的增删改查等操作"""

    def __init__(self):
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
        db_name = "swanlab.sqlite3"
        # 创建数据库文件，如果已经存在则不创建
        db_path = os.path.join(SWANLAB_FOLDER, db_name)
        # 在此处连接数据库，作为一个属性，接下来
        self.__con = sqlite3.connect(db_path)
        # FIXME 如果$expriments数据表单不存在，则创建
        print("swanlab database initialized")

    @property
    def experiments(self) -> List[str]:
        """列出当前数据库中的所有实验"""
        # FIXME 此处应该返回一个列表，包含所有的实验名称
        return []

    def init(
        self,
        expriment_name: str = None,
        description: str = "",
        config: dict = {},
    ):
        """初始化数据库，创建实验级别的数据表单

        Parameters
        ----------
        expriment_name : str, optional
            实验名称，默认自己生成，在内部进行实验名称的排查, by default generate_random_tree_name
        description : str, optional
            实验描述, by default ""
        config : dict, optional
            实验配置, by default {}
        """
        if expriment_name is None:
            expriment_name = generate_random_tree_name()
        # TODO 检查名称是否合法
        check_expriment_name(expriment_name)
        # 获取当前已经存在的实验表单集合
        # 保证实验名称唯一
        expriment_name = make_expriment_name_unique(expriment_name, self.experiments)
        # FIXME 在已经创建的实验表单中添加一条记录

        # FIXME 如果是浮点数，保留4位小数

        return

    def add(self, tag: str, data: Union[str, float], namespace: str = "charts"):
        """添加数据到数据库，保存数据，完成几件事情：
        1. 如果{expriment_name}_{tag}表单不存在，则创建
        2. 添加记录到{expriment_name}_{tag}表单中，包括create_time等
        3. {expriment_name}$chart表单中是否存在此字段，不存在则添加

        Parameters
        ----------
        tag : str
            数据标签，用于区分同一资源下不同的数据
        data : Union[str, float]
            定位到的数据，暂时只支持str和float类型（事实上目前只支持float类型）
        namespace : str, optional
            命名空间，用于区分不同的数据资源（对应{expriment_name}$chart中的tag）, by default "charts"
        """
        pass
