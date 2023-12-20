#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:15:42
@File: swanlab\server\database.py
@IDE: vscode
@Description:
    项目模块，创建项目级别数据库库，接下来针对实验级别的数据在此基础上进行操作
"""
import os
from ..env import swc
from .table import ProjectTablePoxy
from .expriment import ExperimentTable, generate_random_tree_name, check_experiment_name, make_experiment_name_unique
from ..utils import lock_file
from io import TextIOWrapper
import ujson


class ProjectTable(ProjectTablePoxy):
    """实验管理类，用于管理实验，包括创建实验，删除实验，修改实验配置等操作

    # Attributes
    data: dict，实验管理类的数据，json格式
    """

    path = swc.project
    default_data = {"_sum": 0, "experiments": []}

    def __init__(self, data: dict):
        """初始化实验管理类"""
        # print("project init")
        # 保存表单信息
        super().__init__(data, self.path)
        # 当前项目运行的实验
        self.__experiment: ExperimentTable = None

    @property
    def sum(self) -> int:
        """返回实验总数"""
        return self["_sum"]

    @sum.setter
    def sum(self, value):
        """设置实验总数"""
        self["_sum"] = value

    @property
    def experiment(self) -> ExperimentTable:
        """返回当前实验"""
        return self.__experiment

    def add_experiment(self, name: str = None, description: str = None, config: dict = None):
        """添加一个实验

        Parameters
        ----------
        name : str, optional
            实验名称, 如果为None，自动生成一个实验名称
        description : str, optional
            实验描述, 如果为None，设置为`""`
        config : dict, optional
            实验的其他配置，如果为None，设置为`{}`
        """
        # 获取当前已经存在的实验名称集合
        experiments = [item["name"] for item in self["experiments"]]

        # 获取实验名称
        if name is None:
            name = generate_random_tree_name(experiments)
        else:
            check_experiment_name(name)
        # 获取实验描述和配置
        if description is None:
            description = ""
        if config is None:
            config = {}
        # 保证实验名称唯一
        name = make_experiment_name_unique(name, experiments)
        # 创建实验信息
        self.sum = self.sum + 1
        self.__experiment = ExperimentTable(self.sum, name, description, config, len(experiments) + 1)
        # 添加一个实验到self["experiments"]中
        self["experiments"].append(self.__experiment.__dict__())

    @lock_file(file_path=path, mode="r+")
    def success(self, file: TextIOWrapper):
        """实验成功完成，更新实验状态，再次保存实验信息"""
        # 锁上文件，更新实验状态
        project = ujson.load(file)
        self.__experiment.success()
        for index, experiment in enumerate(project["experiments"]):
            if experiment["experiment_id"] == self.__experiment.experiment_id:
                project["experiments"][index] = self.__experiment.__dict__()
                # print("success experiment ", project["experiments"][index])
                break
        self.save(file, project)

    @lock_file(file_path=path, mode="r+")
    def fail(self, file: TextIOWrapper):
        """实验失败，更新实验状态，再次保存实验信息"""
        # 锁上文件，更新实验状态
        project = ujson.load(file)
        self.__experiment.fail()
        for index, experiment in enumerate(project["experiments"]):
            if experiment["experiment_id"] == self.__experiment.experiment_id:
                project["experiments"][index] = self.__experiment.__dict__()
                # print("success experiment ", project["experiments"][index])
                break
        self.save(file, project)
