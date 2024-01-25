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
import datetime
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409, NOT_FOUND_404
from fastapi import Request
from urllib.parse import quote
from ..settings import (
    get_logs_dir,
    get_exp_dir,
    get_config_path,
    get_console_dir,
    get_meta_path,
    get_requirements_path,
)
from ...utils import get_a_lock
from ...utils.file import check_desc_format
from ...utils.time import create_time
from ...utils.font import DEFAULT_COLOR
import yaml
from ...log import swanlog
from typing import List, Dict

from ...db import connect, NotExistedError

from ...db import (
    Project,
    Experiment,
    Chart,
    Namespace,
)

__to_list = Experiment.search2list

# ---------------------------------- 通用 ----------------------------------


# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID
# 实验运行状态
RUNNING_STATUS = Experiment.RUNNING_STATUS


# ---------------------------------- 工具函数 ----------------------------------


def __clear_field(target: List[dict], field: str) -> list[dict]:
    """遍历字典列表清除某个字段

    Parameters
    ----------
    target : List[dict]
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


def __get_exp_dir_by_id(experiment_id: int) -> str:
    """通过 experiment_id 获取实验目录

    Parameters
    ----------
    experiment_id : int
        实验唯一id

    Returns
    -------
    str
        实验目录路径
    """

    return get_exp_dir(Experiment.get(experiment_id).run_id)


def __get_logs_dir_by_id(experiment_id: int) -> str:
    """通过 experiment_id 获取实验 保存目录

    Parameters
    ----------
    experiment_id : int
        实验唯一id

    Returns
    -------
    str
        实验 logs 目录路径
    """
    return get_logs_dir(Experiment.get(experiment_id).run_id)


def __get_console_dir_by_id(experiment_id: int) -> str:
    """通过 experiment_id 获取实验 console 保存目录

    Parameters
    ----------
    experiment_id : int
        实验唯一id

    Returns
    -------
    str
        实验 console 目录路径
    """

    return get_console_dir(Experiment.get(experiment_id).run_id)


def __get_requirements_path_by_id(experiment_id: int):
    """通过 experiment_id 获取实验依赖存储路径

    Parameters
    ----------
    experiment_id : int
        实验唯一id

    Returns
    -------
    str
        实验依赖存储路径
    """

    return get_requirements_path(Experiment.get(experiment_id).run_id)


# ---------------------------------- 路由对应的处理函数 ----------------------------------


# 获取实验信息
def get_experiment_info(experiment_id: int):
    """获取实验信息
    1. 数据库中获取实验的基本信息
    2. 从实验目录获取配置信息
    3. 获取实验元信息

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

    # 加载实验元信息
    meta_path = get_meta_path(experiment["run_id"])
    if os.path.exists(meta_path):
        with get_a_lock(meta_path) as f:
            experiment["system"] = ujson.load(f)

    # 实验默认颜色
    experiment["default_color"] = DEFAULT_COLOR

    return SUCCESS_200(experiment)


# 获取表单数据
def get_tag_data(experiment_id: int, tag: str) -> dict:
    """获取表单数据
    根据实验id得到实验的运行id，然后根据运行id和tag得到实验的数据

    Parameters
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，已经进行了 URIComponent 编码
    """
    # ---------------------------------- 前置处理 ----------------------------------
    # 获取tag对应的存储目录
    try:
        tag_path: str = os.path.join(__get_logs_dir_by_id(experiment_id), tag)
    except NotExistedError:
        return NOT_FOUND_404("experiment not found")
    if not os.path.exists(tag_path):
        return NOT_FOUND_404("tag not found")
    # 获取目录下存储的所有数据
    # 降序排列，最新的数据在最前面
    files: list = os.listdir(tag_path)
    # files中去除_summary.json文件
    files = [f for f in files if not f.endswith("_summary.json")]
    if len(files) == 0:
        return []
    files.sort()
    tag_data: list = []
    # 最后一个文件代表当前数据量
    last_file = files[-1]
    tag_json = None
    # ---------------------------------- 开始读取最后一个文件 ----------------------------------

    # 锁住此文件，不再允许其他进程访问，换句话说，实验日志在log的时候被阻塞
    with get_a_lock(os.path.join(tag_path, last_file), mode="r") as f:
        # 读取数据
        tag_json = ujson.load(f)
        # 倒数第二个文件+当前文件的数据量等于总数据量
        # 倒数第二个文件可能不存在
        count = files[-2].split(".")[0] if len(files) > 1 else 0
        count = int(count) + len(tag_json["data"])
    # 读取完毕，文件解锁
    # ---------------------------------- 返回所有数据 ----------------------------------
    # FIXME: 性能问题
    # 读取所有数据
    # tag_json是最后一个文件的数据
    # 按顺序读取其他文件的数据
    tag_data_list: List[List[Dict]] = []
    for path in files[:-1]:
        # 读取tag数据，由于目前在设计上这些文件不会再被修改，所以不需要加锁
        with open(os.path.join(tag_path, path), "r") as f:
            tag_data_list.append(ujson.load(f)["data"])
    # 将数据合并
    for data in tag_data_list:
        tag_data.extend(data)
    tag_data.extend(tag_json["data"])
    # COMPAT 如果第一个数据没有index，就循环每个数据，加上index
    if tag_data[0].get("index") is None:
        for index, data in enumerate(tag_data):
            data["index"] = str(index + 1)
    # COMPAT 如果第一个数据的index不是string，改为string
    if not isinstance(tag_data[0]["index"], str):
        for data in tag_data:
            data["index"] = str(data["index"])
    # 根据index升序排序
    tag_data.sort(key=lambda x: int(x["index"]))
    # tag_data 的 最后一个数据增加一个字段_last = True
    tag_data[-1]["_last"] = True
    # 获取_summary文件
    summary_path = os.path.join(tag_path, "_summary.json")
    if os.path.exists(summary_path):
        with get_a_lock(summary_path, "r") as f:
            summary = ujson.load(f)
            max_value = summary["max"]
            min_value = summary["min"]
    else:
        # COMPAT 如果_summary文件不存在，手动获取最大值和最小值
        warn = f"Summary file of tag '{tag}' not found, SwanLab will automatically get the maximum and minimum values."
        swanlog.warning(warn)
        # 遍历tag_data，获取最大值和最小值
        # 提取 data 字段的值
        data_values = [entry["data"] for entry in tag_data]
        # 获取最大值和最小值
        max_value = max(data_values)
        min_value = min(data_values)
    return SUCCESS_200(data={"sum": len(tag_data), "max": max_value, "min": min_value, "list": tag_data})


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
    for index, chart in enumerate(charts):
        sources = []
        for source in __to_list(chart.sources):
            sources.append(source["tag_id"]["name"])
        chart_list[index]["source"] = sources

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
        summaries: list
            每个tag的最后一个数据
    """

    experiment = Experiment.get_by_id(experiment_id)
    tag_list = [tag["name"] for tag in __to_list(experiment.tags)]
    experiment_path = __get_logs_dir_by_id(experiment_id)
    tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
    summaries = []
    for tag in tag_list:
        if quote(tag, safe="") not in tags:
            summaries.append({"key": tag, "value": "TypeError"})
            continue
        tag_path = os.path.join(experiment_path, quote(tag, safe=""))
        logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
        with get_a_lock(os.path.join(tag_path, logs[-1]), mode="r") as f:
            data = ujson.load(f)
            # str 转化的目的是为了防止有些不合规范的数据导致返回体对象化失败
            data = str(data["data"][-1]["data"])
            summaries.append({"key": tag, "value": data})

    return SUCCESS_200({"summaries": summaries})


MAX_NUM = 6000


# 获取实验最近日志
def get_recent_logs(experiment_id):
    """一下返回最多 MAX_NUM 条打印记录

    Parameters
    ----------
    experiment_id : int
        实验唯一ID

    Returns
    -------
    dict :
        recent: List[str]
            0: 截止处日志文件
            1: 开始处日志文件
        logs: List[str]
            日志列表
    """

    console_path: str = __get_console_dir_by_id(experiment_id)
    consoles: list = [f for f in os.listdir(console_path)]
    # 含有error.log，在返回值中带上其中的错误信息
    error = None
    if "error.log" in consoles:
        with open(os.path.join(console_path, "error.log"), mode="r") as f:
            error = f.read().split("\n")
        # 在consoles里删除error.log
        consoles.remove("error.log")
    total: int = len(consoles)
    # # 如果 total 大于 1, 按照时间排序
    if total > 1:
        consoles = sorted(consoles, key=lambda x: datetime.strptime(x[:-4], "%Y-%m-%d"), reverse=True)
    logs = []
    # # current_page = total
    for index, f in enumerate(consoles, start=1):
        with open(os.path.join(console_path, f), mode="r", encoding="utf-8") as log:
            logs.extend(log.read().split("\n"))
            # 如果当前收集到的数据超过限制，退出循环
            if len(logs) >= MAX_NUM:
                # current_page = index
                break
    logs = logs[:MAX_NUM]
    end = (logs[-1] if not logs[-1] == "" else logs[-2]).split(" ")[0]
    data = {
        "recent": [logs[0].split(" ")[0], end],
        "logs": logs,
    }
    if error is not None:
        data["error"] = error
    # 返回最新的 MAX_NUM 条记录
    return SUCCESS_200(data)


# 获取图表信息
def get_experimet_charts(experiment_id: int):
    """获取图表信息

    Parameters
    ----------
    experiment_id : int
        实验唯一 ID

    Returns
    -------
    dict :
        _sum: integer
        charts: List[dict]
        namesapces: List[dict]
    """

    charts = Chart.filter(Chart.experiment_id == experiment_id)
    chart_list = __to_list(charts)
    # 获取每个图表对应的数据源
    for index, chart in enumerate(charts):
        sources = []
        for source in __to_list(chart.sources):
            sources.append(source["tag_id"]["name"])
            if source["error"] is not None and source["error"] != "":
                try:
                    chart_list[index]["error"] = ujson.loads(source["error"])
                except:
                    pass
        chart_list[index]["source"] = sources

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
            "_sum": charts.count(),
            "charts": __clear_field(chart_list, "experiment_id"),
            "namespaces": namespace_list,
        }
    )


# 更改项目信息
async def update_experiment_info(experiment_id: int, request: Request):
    """修改实验的信息

    Parameters
    ----------
    experiment_id : int
        实验id
    body : Body
        name: str
            实验名称
        description: str
            实验描述

    Returns
    -------
    dict :
        name: str
        description: str
    """

    db = connect()
    body = await request.json()
    with db.atomic():
        experiment = Experiment.get(experiment_id)
        experiment.name = body.get("name")
        experiment.description = body.get("description")
        experiment.save()

    db.close()

    return SUCCESS_200(body)


# 删除实验
def delete_experiment(experiment_id: int):
    """删除实验

    注意，需要先判断当前实验是否正在运行中，不可删除运行中的实验
    同时，需要删除所有表中含有该 experiment_id 的行

    Parameters
    ----------
    experiment_id : Int
        实验唯一ID

    Returns
    -------
    project : Dictionary
        删除实验后的项目信息，提供给前端更新界面

    TODO: 该怎么删才能高效删干净
    """

    # 先删除实验目录
    experiment_path = __get_exp_dir_by_id(experiment_id)
    shutil.rmtree(experiment_path)
    # 再清除数据库中的实验数据
    Experiment.delete().where(Experiment.id == experiment_id).execute()

    return SUCCESS_200({"experiment_id": experiment_id})


# 停止实验
def stop_experiment(experiment_id: int):
    """停止实验

    Parameters
    ----------
    experiment_id : Int
        实验唯一ID

    Returns
    -------
    project : Dictionary
        停止实验后的项目信息，提供给前端更新界面
    """

    db = connect()
    with db.atomic():
        experiment = Experiment.get(experiment_id)
        experiment.status = Experiment.STOPPED_STATUS
        experiment.save()

    return SUCCESS_200(
        {
            "id": experiment_id,
            "status": Experiment.STOPPED_STATUS,
            "update_time": create_time(),
        }
    )


# 获取实验依赖
def get_experiment_requirements(experiment_id: int):
    """获取实验依赖

    Parameters
    ----------
    experiment_id : int
        实验唯一ID

    Returns
    -------
    dict :
        requirements: list
            每个依赖项为一行，以列表的形式返回
    """

    path = __get_requirements_path_by_id(experiment_id)
    if not os.path.exists(path):
        return DATA_ERROR_500("failed to find requirements")
    with open(path) as f:
        requirements = f.read()

    return SUCCESS_200({"requirements": requirements.split("\n")})
