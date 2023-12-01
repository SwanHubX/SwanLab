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
import os
from ..env import SWANLAB_LOGS_FOLDER


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
        self.index = index

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
            "description": self.description,
            "config": self.config,
            "index": self.index,
            "tags": self.tags,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }
