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
import json
import os
from ..env import SWANLAB_FOLDER


class ExperimentTable(ExperimentPoxy):
    """初始化一个实验，需要传入一些实验的基本信息，如实验名称，实验描述，实验配置等信息"""

    def __init__(self, name: str, description: str, config: dict, sort: int):
        # 初始化一个实验配置，name必须保证唯一，但是不在此处检查，而是在创建实验的时候检查
        super().__init__(name, path=os.path.join(SWANLAB_FOLDER, "logs"))
        self.name = name
        self.tags = []
        self.description = description
        self.config = config
        self.sort = sort

    def __dict__(self) -> dict:
        """序列化此对象

        Returns
        -------
        dict
            序列化后的字典对象
        """
        return {
            "name": self.name,
            "description": self.description,
            "config": self.config,
            "sort": self.sort,
            "tags": self.tags,
            "create_time": self.create_time,
            "update_time": self.update_time,
        }
