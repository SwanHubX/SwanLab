#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:15:42
@File: swanlab\server\database.py
@IDE: vscode
@Description:
    项目模块，创建项目级别数据库库，接下来针对实验级别的数据在此基础上进行操作
"""
from typing import Union
import os
from ..env import SWANLAB_FOLDER
from .experiments_name import generate_random_tree_name, check_experiment_name, make_experiment_name_unique
import ujson
from datetime import datetime
from .table import ProjectTablePoxy
from .expriment import ExperimentTable

SWANLAB_LOGS_FOLDER = os.path.join(SWANLAB_FOLDER, "logs")

DEFAULT_CHART = {
    "chart_id": 0,
    "tag": "default",
    "source": [],
    "type": "default",
    "config": {},
    "created_time": datetime.now().isoformat(),
}


class ProjectTable(ProjectTablePoxy):
    """实验管理类，用于管理实验，包括创建实验，删除实验，修改实验配置等操作

    # Attributes
    data: dict，实验管理类的数据，json格式
    """

    path = os.path.join(SWANLAB_LOGS_FOLDER, "project.json")

    def __init__(self):
        """初始化实验管理类"""
        project_exist = os.path.exists(self.path)
        # 判断path是否存在，如果存在，则加载数据，否则创建
        if project_exist:
            with open(self.path, "r") as f:
                data = ujson.load(f)
        else:
            data = {"_sum": 0, "experiments": []}
        # 保存表单信息
        super().__init__(data, self.path)
        # 保存表单信息
        if not project_exist:
            self.save()
        self._experiment: ExperimentTable = None

    @property
    def sum(self):
        """返回实验总数"""
        return self["_sum"]

    @property
    def expriments(self):
        """列出当前数据库中的所有实验"""
        experiments = self["experiments"]
        # 此处应该返回一个列表，包含所有的实验名称
        return [item["name"] for item in experiments]

    @property
    def add_expriment(self, name: str = None, description: str = None, config: dict = None):
        """添加一个实验"""
        if name is None:
            name = generate_random_tree_name()
        else:
            check_experiment_name(name)
        if description is None:
            description = ""
        if config is None:
            config = {}
        # 获取当前已经存在的实验名称集合
        experiments = [item["name"] for item in self.data["experiments"]]
        # 保证实验名称唯一
        name = make_experiment_name_unique(name, experiments)
        # 创建实验信息
        self._experiment = ExperimentTable(name, description, config, self.sum)


class ChartTable(ProjectTablePoxy):
    """图表管理类，用于管理图表，包括创建图表，删除图表，修改图表配置等操作"""

    path = os.path.join(SWANLAB_LOGS_FOLDER, "chart.json")

    def __init__(self):
        """初始化图表管理类"""
        # 判断path是否存在，如果存在，则加载数据，否则创建
        if os.path.exists(self.path):
            with open(self.path, "r") as f:
                data = ujson.load(f)
        else:
            data = {
                "charts": [],
            }
        # 保存表单信息
        super().__init__(data, self.path)
        # 保存表单信息
        self.save()


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
        # 此时必须保证.swanlab文件夹存在，但是这并不是本类的指责，所以不检查
        # 需要检查logs文件夹是否存在，不存在则创建
        if not os.path.exists(SWANLAB_LOGS_FOLDER):
            os.mkdir(SWANLAB_LOGS_FOLDER)
        # 如果ProjectTable.path都不存在或者都存在则不报错
        project_exist = os.path.exists(ProjectTable.path)
        chart_exist = os.path.exists(ChartTable.path)
        if not project_exist and not chart_exist or project_exist and chart_exist:
            pass
        else:
            # 如果只存在一个，则报错，这是不可能的
            raise FileExistsError("invalid project logs")
        # 项目基础表单
        self._project: ProjectTable = ProjectTable()
        # 图表基础表单
        self._chart: ChartTable = ChartTable()

    def init(
        self,
        experiment_name: str = None,
        description: str = "",
        config: dict = {},
    ):
        """初始化数据库，创建实验级别的数据表单

        Parameters
        ----------
        experiment_name : str, optional
            实验名称，默认自己生成，在内部进行实验名称的排查, by default generate_random_tree_name
        description : str, optional
            实验描述, by default ""
        config : dict, optional
            实验配置, by default {}
        """
        # 获取之前的实验记录
        with open(expriments_path, "r") as f:
            expriments_json = ujson.load(f)
        # 写入新的实验
        with open(expriments_path, "w") as f:
            expriments_json["__index"] += 1
            expriments_json["experiments"].append(
                {
                    "expriment_id": expriments_json["__index"],
                    "name": experiment_name,
                    "index": 0,
                    "description": description,
                    "config": config,
                    "create_time": datetime.now().isoformat(),
                }
            )
            # 保存实验id
            self.experiment_id = expriments_json["experiments"][-1]["expriment_id"]
            ujson.dump(expriments_json, f)
        # 创建实验专属目录
        if not os.path.exists(os.path.join(SWANLAB_FOLDER, experiment_name)):
            os.makedirs(os.path.join(SWANLAB_FOLDER, experiment_name))
        # 实验表格配置
        chart_json = {
            "__index": 0,
            "charts": [DEFAULT_CHART],
            "create_time": datetime.now().isoformat(),
            "update_time": datetime.now().isoformat(),
            "experiment_id": self.experiment_id,
        }
        chart_path = os.path.join(SWANLAB_FOLDER, experiment_name, EXPERIMENT_CHART)
        if os.path.exists(chart_path):
            with open(chart_path, "r") as f:
                chart_json = ujson.load(f)
        else:
            with open(chart_path, "w") as f:
                ujson.dump(chart_json, f)
        # FIXME 如果是浮点数，保留4位小数
        return

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
