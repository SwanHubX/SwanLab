#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-01 19:33:40
@File: swanlab/database/expriment.py
@IDE: vscode
@Description:
    实验类，操作和记录实验相关的信息
"""
from .table import ExperimentPoxy
from .chart import ChartTable
import os
from ..env import SWANLAB_LOGS_FOLDER
from typing import Union
from ..utils import create_time
import sys


class ExperimentTable(ExperimentPoxy):
    """初始化一个实验，需要传入一些实验的基本信息，如实验名称，实验描述，实验配置等信息"""

    def __init__(self, experiment_id: int, name: str, description: str, config: dict, index: int):
        # 初始化一个实验配置，name必须保证唯一，但是不在此处检查，而是在创建实验的时候检查
        # 创建name对应的文件夹
        path = os.path.join(SWANLAB_LOGS_FOLDER, name)
        if not os.path.exists(path):
            os.mkdir(path)
        super().__init__(path)
        self.experiment_id = experiment_id
        self.name = name
        self.tags = []
        self.description = description
        self.config = config
        self.argv = sys.argv
        self.index = index
        self.status = 0  # 0: 正在运行，1: 运行成功，-1: 运行失败
        self.__chart = ChartTable(base_path=path)

    def __dict__(self) -> dict:
        """序列化此对象

        Returns
        -------
        dict
            序列化后的字典对象
        """
        return {
            "experiment_id": self.experiment_id,
            "name": self.name,
            "status": self.status,
            "description": self.description,
            "config": self.config,
            "argv": self.argv,
            "index": self.index,
            "tags": self.tags,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }

    def is_tag_exist(self, tag: str) -> bool:
        """判断tag是否已经存在于此实验中

        Parameters
        ----------
        tag : str
            tag名称

        Returns
        -------
        bool
            如果存在，返回True，否则返回False
        """
        return tag in self.tags

    def add(self, tag: str, data: Union[float, str], namespace: str):
        """记录一条新的tag数据

        Parameters
        ----------
        tag : str
            tag名称
        data : Union[float, str]
            tag数据，可以是浮点型，也可以是字符串，但不可以是其他类型
        namespace : str
            命名空间，用于区分不同的数据资源（对应{experiment_name}$chart中的tag）
        """
        if not self.is_tag_exist(tag):
            # 在chart中记录此tag
            # self.__chart.add(tag=tag)
            # 在实验中记录此tag
            self.tags.append(tag)
        # 往本地添加新的数据
        self.save_tag(tag, data, self.experiment_id)

    def success(self):
        """实验成功完成，更新实验状态"""
        self.update_time = create_time()
        self.status = 1
