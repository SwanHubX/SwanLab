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
from ..env import SWANLAB_FOLDER
from .experiments_name import generate_random_tree_name, check_experiment_name, make_experiment_name_unique
import ujson
from datetime import datetime

PROJECT_CONFIG = "experiments.json"
EXPERIMENT_CHART = "chart.json"
DEFAULT_CONFIG = {
    "__index": 0,
    "experiments": [],
    "create_time": datetime.now().isoformat(),
    "update_time": datetime.now().isoformat(),
}
DEFAULT_CHART = {
    "chart_id": 0,
    "tag": "default",
    "source": [],
    "type": "default",
    "config": {},
    "created_time": datetime.now().isoformat(),
}


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
        if not os.path.exists(SWANLAB_FOLDER):
            os.makedirs(SWANLAB_FOLDER)
        # 检查项目的基础配置
        expriments_path = os.path.join(SWANLAB_FOLDER, PROJECT_CONFIG)
        if not os.path.exists(expriments_path) or os.path.getsize(expriments_path) == 0:
            # 创建 experiments.json 文件
            with open(expriments_path, "w") as f:
                ujson.dump(DEFAULT_CONFIG, f)

        print("swanlab database initialized")

    @property
    def experiments(self) -> List[str]:
        """列出当前数据库中的所有实验"""
        experiments = []
        # 此处应该返回一个列表，包含所有的实验名称
        with open(os.path.join(SWANLAB_FOLDER, PROJECT_CONFIG)) as f:
            experiments = ujson.load(f)["experiments"]
        return [item["name"] for item in experiments]

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
        if experiment_name is None:
            experiment_name = generate_random_tree_name()
        # TODO 检查名称是否合法
        check_experiment_name(experiment_name)
        # 获取当前已经存在的实验表单集合
        # 保证实验名称唯一
        experiment_name = make_experiment_name_unique(experiment_name, self.experiments)
        # 给实例绑定实验名
        self.experiment_name = experiment_name

        # 创建实验目录及实验配置
        expriments_json = {}
        expriments_path = os.path.join(SWANLAB_FOLDER, PROJECT_CONFIG)
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
        # 如果json不存在或者为空
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
