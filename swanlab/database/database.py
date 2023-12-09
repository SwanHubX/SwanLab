#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-02 00:38:37
@File: swanlab\database\database.py
@IDE: vscode
@Description:
    数据库模块，连接project表单对象
"""
import os
from ..env import SWANLAB_LOGS_FOLDER
from .project import ProjectTable
from .expriment import ExperimentTable
from ..utils import lock_file
from typing import Union
from io import TextIOWrapper
import ujson
import atexit

# flag，代表已经执行了inited函数
inited = False


class SwanDatabase(object):
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
        # 此时必须保证.swanlab文件夹存在，但是这并不是本类的职责，所以不检查

        # TODO 但是目前还是先在这里创建
        from ..env import SWANLAB_FOLDER

        if not os.path.exists(SWANLAB_FOLDER):
            os.mkdir(SWANLAB_FOLDER)

        # 需要检查logs文件夹是否存在，不存在则创建
        if not os.path.exists(SWANLAB_LOGS_FOLDER):
            os.mkdir(SWANLAB_LOGS_FOLDER)
        # 项目基础表单
        self.__project: ProjectTable = None
        # 如果项目配置文件不存在，创建
        open(ProjectTable.path, "a", encoding="utf-8").close()
        # 表单会在init中创建，所有的创建会在一个文件读取周期内完成，以防止多进程写入同一个文件带来的问题

    @lock_file(file_path=ProjectTable.path, mode="r+")
    def init(
        self,
        file: TextIOWrapper,
        experiment_name: str = None,
        description: str = "",
        config: dict = {},
    ):
        """初始化项目级别表单，创建实验级别的数据表单

        Parameters
        ----------
        experiment_name : str, optional
            实验名称，默认自己生成，在内部进行实验名称的排查, by default generate_random_tree_name
        description : str, optional
            实验描述, by default ""
        config : dict, optional
            实验配置, by default {}
        file : TextIOWrapper, optional
            文件对象，用于文件锁定, by default None
        """
        # 同一运行时不允许运行两次此函数，通过flag来实现
        global inited
        if inited:
            raise RuntimeError("Swanlab has already been inited!")
        inited = True

        # 创建回调函数函数
        def callback():
            sd.success()

        # 注册此函数
        atexit.register(callback)
        # 检查实验名称是否存在
        project_exist = os.path.exists(ProjectTable.path) and os.path.getsize(ProjectTable.path) != 0
        # 初始化项目对象
        self.__project = ProjectTable(data=ujson.load(file) if project_exist else ProjectTable.default_data)
        # 创建实验
        self.__project.add_experiment(experiment_name, description, config)
        # 保存项目表单
        self.__project.save(file)
        # 在修饰器中自动解锁

    @property
    def experiment(self) -> ExperimentTable:
        """获取当前实验对象"""
        return self.__project.experiment

    def add(self, tag: str, data: Union[str, float], namespace: str = "charts"):
        """添加数据到数据库，保存数据，完成几件事情：
        1. 如果{experiment_name}_{tag}表单不存在，则创建
        2. 添加记录到{experiment_name}_{tag}表单中，包括create_time等
        3. chart表单中是否存在此字段，不存在则添加
        当然这些工作并不由本函数完成，而是由表单对象完成

        Parameters
        ----------
        tag : str
            数据标签，用于区分同一资源下不同的数据
        data : Union[str, float]
            定位到的数据，暂时只支持str和float类型（事实上目前只支持float类型）
        namespace : str, optional
            命名空间，用于区分不同的数据资源（对应{experiment_name}$chart中的tag）, by default "charts"
        """
        global inited
        if not inited:
            raise RuntimeError("Swanlab must be inited first!")
        self.__project.experiment.add(tag, data, namespace)

    def success(self):
        """标记实验成功"""
        self.__project.success()


sd = SwanDatabase()
