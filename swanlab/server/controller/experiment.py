#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-20 21:21:45
@File: swanlab\server\controller\experiment.py
@IDE: vscode
@Description:
    实验相关 api 的处理函数
"""

import os
import ujson
import shutil
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409
from fastapi import Request
from urllib.parse import unquote
from ..settings import (
    get_logs_dir,
    get_tag_dir,
    get_files_dir,
    get_exp_dir,
    DB_PATH,
    get_config_path,
)
from ...utils import get_a_lock
from ...utils.file import check_desc_format
import yaml

from ...db import (
    tables,
    connect,
)

from ...db import (
    Project,
    Experiment,
    Tag,
    Chart,
    Namespace,
    Display,
)

__to_list = Experiment.search2list

# ---------------------------------- 通用 ----------------------------------


# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID
# 实验运行状态
RUNNING_STATUS = Experiment.RUNNING_STATUS


# ---------------------------------- 工具函数 ----------------------------------


def __clear_field(target: list[dict], field: str) -> list[dict]:
    """遍历字典列表清除某个字段

    Parameters
    ----------
    target : list[dict]
        需要处理的列表
    field : str
        需要删除的字段

    Returns
    -------
    list[dict]
        处理后的字典列表
    """

    for item in target:
        item.pop(field)

    return target


def __get_logs_dir_by_id(experiment_id: int) -> str:
    """通过 experiment_id 获取实验 保存目录

    Parameters
    ----------
    experiment_id : int
        实验唯一id

    Returns
    -------
    str
        实验 run_id
    """

    return get_logs_dir(Experiment.get(experiment_id).run_id)


# ---------------------------------- 路由对应的处理函数 ----------------------------------


# 获取实验信息
def get_experiment_info(experiment_id: int):
    """获取实验信息
    1. 数据库中获取实验的基本信息
    2. 从实验目录获取配置信息

    Parameters
    ----------
    experiment_id : int
        实验唯一id
    """

    experiment = Experiment.get(experiment_id).__dict__()
    experiment.pop("project_id")

    # 加载实验配置
    config_path = get_config_path(experiment["run_id"])
    if os.path.exists(config_path):
        with get_a_lock(config_path) as f:
            experiment["config"] = yaml.load(f, Loader=yaml.FullLoader)

    return SUCCESS_200({"experiment": experiment})


# 获取表单数据
def get_tag_data(experiment_id: int, tag: str) -> dict:
    """获取表单数据
    根据实验id得到实验的运行id，然后根据运行id和tag得到实验的数据
    """

    return SUCCESS_200({})


# 获取实验状态
def get_experiment_status(experiment_id: int):
    """获取实验状态以及实验图表配置，用于实时更新实验状态

    Parameters
    ----------
    experiment_id : int
        实验唯一id，路径传参
    """

    experiment = Experiment.get(experiment_id)
    charts = Chart.filter(Chart.experiment_id == experiment_id)
    chart_list = __clear_field(__to_list(charts), "experiment_id")

    # 当前实验下的命名空间
    namespaces = Namespace.filter(Namespace.experiment_id == experiment_id)
    namespace_list = __clear_field(__to_list(namespaces), "experiment_id")
    # 获取每个命名空间对应的 display
    # display 含有 chart 与 namespace 的对应关系
    for index, namespace in enumerate(namespaces):
        displays = []
        for display in __to_list(namespace.displays):
            displays.append(display["chart_id"]["id"])
        namespace_list[index]["charts"] = displays

    return SUCCESS_200(
        {
            "status": experiment.status,
            "update_time": experiment.update_time,
            "charts": {
                "_sum": charts.count(),
                "charts": chart_list,
                "namespaces": namespace_list,
            },
        }
    )


# 获取实验的总结数据
def get_experiment_summary(experiment_id: int) -> dict:
    """获取实验的总结数据——每个tag的最后一个setp的data

    Parameters
    ----------
    experiment_id : int
        实验id

    Returns
    -------
    dict:
        summaries: array
            每个tag的最后一个数据
    """

    experiment_path = __get_logs_dir_by_id(experiment_id)
    tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
    summaries = []
    for tag in tags:
        tag_path = os.path.join(experiment_path, tag)
        logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
        with get_a_lock(os.path.join(tag_path, logs[-1]), mode="r") as f:
            data = ujson.load(f)
            data = data["data"][-1]["data"]
            summaries.append({"key": unquote(tag), "value": data})
    return SUCCESS_200({"summaries": summaries})
