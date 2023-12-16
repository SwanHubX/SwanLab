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
from ..env import swc
from typing import Union
from ..utils import create_time, generate_color
from .system import get_system_info
import sys


class ExperimentTable(ExperimentPoxy):
    """初始化一个实验，需要传入一些实验的基本信息，如实验名称，实验描述，实验配置等信息"""

    def __init__(self, experiment_id: int, name: str, description: str, config: dict, index: int):
        # 初始化一个实验配置，name必须保证唯一，但是不在此处检查，而是在创建实验的时候检查
        # 创建name对应的文件夹
        swc.add_exp(name)
        if not os.path.exists(swc.logs_folder):
            os.makedirs(swc.logs_folder)
        super().__init__(swc.logs_folder)
        self.experiment_id = experiment_id
        self.name = name
        # tags数据不会被序列化
        self.tags = []
        self.system = get_system_info()
        self.description = description
        self.config = config
        self.argv = sys.argv
        self.index = index
        self.status = 0  # 0: 正在运行，1: 运行成功，-1: 运行失败
        self.__chart = ChartTable(experiment_id=experiment_id)
        self.color = generate_color(experiment_id)

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
            "color": self.color,
            "system": self.system,
            "argv": self.argv,
            "index": self.index,
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
        return any(item["tag"] == tag for item in self.tags)

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
            self.__chart.add(tag=tag)
            # 在实验中记录此tag
            self.tags.append({"tag": tag, "num": 0})
        # 更新tag的数量，并拿到tag的索引
        tag_index = self.update_tag_num(tag)
        # 往本地添加新的数据
        self.save_tag(tag, data, self.experiment_id, tag_index)

    def update_tag_num(self, tag: str) -> int:
        for index, item in enumerate(self.tags):
            if item["tag"] == tag:
                item["num"] += 1
                return item["num"]
        raise ValueError(f"tag={tag} not exist in experiment={self.name}")

    def success(self):
        """实验成功完成，更新实验状态"""
        self.update_time = create_time()
        self.status = 1

    def fail(self):
        """实验失败，更新状态"""
        self.update_time = create_time()
        self.status = -1
