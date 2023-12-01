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
from ..env import SWANLAB_LOGS_FOLDER, SWANLAB_FOLDER
from .project import ProjectTable
from .chart import ChartTable
import portalocker
from typing import Union
import ujson


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
        # 此时必须保证.swanlab文件夹存在，但是这并不是本类的指责，所以不检查
        # 需要检查logs文件夹是否存在，不存在则创建
        if not os.path.exists(SWANLAB_LOGS_FOLDER):
            os.mkdir(SWANLAB_LOGS_FOLDER)
        # 项目基础表单
        self._project: ProjectTable = None
        # 图表基础表单
        self._chart: ChartTable = None
        # 表单会在init中创建，所有的创建会在一个文件读取周期内完成，以防止多进程写入同一个文件带来的问题

    def init(
        self,
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
        """
        # 检查实验名称是否存在
        project_exist = os.path.exists(ProjectTable.path) and os.path.getsize(ProjectTable.path) != 0
        # 项目表单
        f = open(ProjectTable.path, "r+") if project_exist else open(ProjectTable.path, "w+")
        # 锁定文件
        portalocker.lock(f, portalocker.LOCK_EX)
        # 初始化项目对象
        self._project = ProjectTable(data=ujson.load(f) if project_exist else ProjectTable.default_data)
        # 创建实验
        self._project.add_experiment(experiment_name, description, config)
        # 保存项目表单，释放文件锁
        self._project.save(f)
        f.close()

    def add(self, tag: str, data: Union[str, float], namespace: str = "charts"):
        """添加数据到数据库，保存数据，完成几件事情：
        1. 如果{experiment_name}_{tag}表单不存在，则创建
        2. 添加记录到{experiment_name}_{tag}表单中，包括create_time等
        3. {experiment_name}$chart表单中是否存在此字段，不存在则添加

        Parameters
        ----------
        tag : str
            数据标签，用于区分同一资源下不同的数据
        data : Union[str, float]
            定位到的数据，暂时只支持str和float类型（事实上目前只支持float类型）
        namespace : str, optional
            命名空间，用于区分不同的数据资源（对应{experiment_name}$chart中的tag）, by default "charts"
        """
        database = os.path.join(SWANLAB_FOLDER, self.experiment_name, f"{self.experiment_name}_{tag}.json")
        previous_data = {}
        # 如果对应的database文件不存在，则创建
        if not os.path.exists(database) or os.path.getsize(database) == 0:
            previous_data = {
                "data": [],
                "create_time": datetime.now().isoformat(),
                "update_time": datetime.now().isoformat(),
            }
        else:
            with open(database, "r") as f:
                previous_data = ujson.load(f)
        # 添加记录
        previous_data["data"].append(
            {
                "create_time": datetime.now().isoformat(),
                "data": data,
                "tag_id": len(previous_data["data"]),
                "expriment_id": self.experiment_id,
            }
        )
        with open(database, "w") as f:
            ujson.dump(previous_data, f)
