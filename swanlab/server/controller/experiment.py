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
from datetime import datetime
from ..module.resp import SUCCESS_200, DATA_ERROR_500, NOT_FOUND_404
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
from ...utils.time import create_time
import yaml
from ...log import swanlog
from .db import (
    Project,
    Experiment,
    Chart,
    Tag,
    connect,
    NotExistedError,
)
from .utils import get_exp_charts, clear_field, read_tag_data, get_tag_files, LOGS_CONFIGS, lttb
from ...compat.server.controller.experiment import compat_text

__to_list = Experiment.search2list

# ---------------------------------- 通用 ----------------------------------


# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID
# 实验运行状态
RUNNING_STATUS = Experiment.RUNNING_STATUS


# ---------------------------------- 工具函数 ----------------------------------


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
        with open(config_path, 'r+') as f:
            experiment["config"] = yaml.load(f, Loader=yaml.FullLoader)

    # 加载实验元信息
    meta_path = get_meta_path(experiment["run_id"])
    if os.path.exists(meta_path) and not os.stat(meta_path).st_size == 0:
        with open(meta_path, 'r+') as f:
            experiment["system"] = ujson.load(f)
    else:
        experiment["system"] = {}
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
    tag_data: list = []
    # ---------------------------------- 读取文件数据 ----------------------------------
    current_logs = get_tag_files(tag_path, LOGS_CONFIGS)
    for file in current_logs:
        tag_data = tag_data + read_tag_data(os.path.join(tag_path, file))
    # ---------------------------------- 返回所有数据 ----------------------------------
    # 如果数据为空，返回空列表
    if len(tag_data) == 0:
        return SUCCESS_200(data={"sum": 0, "list": [], "experiment_id": experiment_id})
    # COMPAT 如果第一个数据没有index，就循环每个数据，加上index
    if tag_data[0].get("index") is None:
        for index, data in enumerate(tag_data):
            data["index"] = str(index + 1)
    # COMPAT 如果第一个数据的index不是int，改为int
    if not isinstance(tag_data[0]["index"], int):
        for data in tag_data:
            data["index"] = int(data["index"])
    # COMPAT 如果当前tag的类型为text，并且在media文件夹下存在相同名文件夹
    texts = compat_text(experiment_id, tag)
    # 如果需要兼容，进行重新赋值，直接获取文本内容
    if texts:
        tag_data = texts
    # 根据index升序排序
    tag_data.sort(key=lambda x: int(x["index"]))
    # tag_data 的 最后一个数据增加一个字段_last = True
    tag_data[-1]["_last"] = True
    # 获取_summary文件
    summary_path = os.path.join(tag_path, "_summary.json")
    if os.path.exists(summary_path):
        with open(summary_path, "r") as f:
            summary = ujson.load(f)
            max_value = summary.get("max", None)
            min_value = summary.get("min", None)
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
    return SUCCESS_200(
        data={
            "sum": len(tag_data),
            "max": max_value,
            "min": min_value,
            "list": lttb(tag_data),
            # 标注此数据隶属于哪个实验
            "experiment_id": experiment_id,
        }
    )


# 获取实验状态
def get_experiment_status(experiment_id: int):
    """获取实验状态以及实验图表配置，用于实时更新实验状态

    Parameters
    ----------
    experiment_id : int
        实验唯一id，路径传参
    """

    experiment = Experiment.get(experiment_id)
    chart_list, namespace_list = get_exp_charts(experiment_id)

    return SUCCESS_200(
        {
            "status": experiment.status,
            "update_time": experiment.update_time,
            "finish_time": experiment.finish_time,
            "charts": {
                "charts": clear_field(chart_list, "experiment_id"),
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
    # 通过外键反链获取实验下的所有tag
    tag_list = [tag["name"] for tag in __to_list(experiment.tags)]
    experiment_path = __get_logs_dir_by_id(experiment_id)
    # 通过目录结构获取所有正常的tag
    tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
    # 实验总结数据
    summaries = []
    for tag in tag_list:
        # 如果 tag 记录在数据库，但是没有对应目录，说明 tag 有问题
        # 所以 tags 是 tag_list 的子集，出现异常的 tag 会记录在数据库但不会添加到目录结构中
        if quote(tag, safe="") not in tags:
            summaries.append({"key": tag, "value": "TypeError"})
            continue
        tag_path = os.path.join(experiment_path, quote(tag, safe=""))
        # 获取 tag 目录下的所有存储的日志文件
        logs = get_tag_files(tag_path, LOGS_CONFIGS)
        # 打开 tag 目录下最后一个存储文件，获取最后一条数据
        with open(os.path.join(tag_path, logs[-1])) as f:
            lines = f.readlines()
            # 最后一行数据，如果为空，取倒数第二行
            data = lines[-1] if len(lines[-1]) else lines[-2]
            last_data = ujson.loads(data)
            summaries.append({"key": tag, "value": last_data["data"]})
    # 获取数据库记录时在实验下的排序
    sorts = {item["name"]: item["sort"] for item in __to_list(experiment.tags)}
    # COMPAT 如果 sorts 中的值不都为 0，说明是新版添加排序后的 tag，这里进行排序 (如果是旧版没有排序的tag，直接按照数据库顺序即可)
    if not all(value == 0 for value in sorts.values()):
        temp = [0] * len(summaries)
        for item in summaries:
            temp[sorts[item["key"]]] = item
        summaries = temp
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
    # 是否存在
    if not os.path.exists(console_path):
        return NOT_FOUND_404("Log Folder not Found")
    consoles: list = [f for f in os.listdir(console_path)]
    # 含有error.log，在返回值中带上其中的错误信息
    error = None
    if "error.log" in consoles:
        with open(os.path.join(console_path, "error.log"), mode="r") as f:
            error = f.read().split("\n")
        # 在consoles里删除error.log
        consoles.remove("error.log")
    # 没有日志，也没有错误日志
    if len(consoles) == 0 and error is None:
        return NOT_FOUND_404("No Logs Found")
    # 如果只有错误日志，直接返回
    elif len(consoles) == 0 and error:
        return SUCCESS_200({"error": error})

    total: int = len(consoles)
    # # 如果 total 大于 1, 按照时间排序
    if total > 1:
        consoles = sorted(consoles, key=lambda x: datetime.strptime(x[:-4], "%Y-%m-%d"), reverse=True)
    logs = []
    # # current_page = total
    for index, f in enumerate(consoles, start=1):
        with open(os.path.join(console_path, f), mode="r", encoding="utf-8") as log:
            logs = log.read().split("\n") + logs
            # 如果当前收集到的数据超过限制，退出循环
            # if len(logs) >= MAX_NUM:
            #     # current_page = index
            #     break
    # 如果 logs 内容为空
    if len(logs) == 0:
        return NOT_FOUND_404("No Logs Found")

    # logs = logs[:MAX_NUM]
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
    chart_list, namespace_list = get_exp_charts(experiment_id)

    return SUCCESS_200(
        {
            "charts": clear_field(chart_list, "experiment_id"),
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
    body["name"] = body["name"].strip()
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
    # 看看该实验下有哪儿些 tag
    tags = Tag.filter(Tag.experiment_id == experiment_id)
    tag_names = [tag["name"] for tag in __to_list(tags)]
    # 检查多实验图表是否有需要删除的
    project_id = Experiment.get_by_id(experiment_id).project_id.id
    # 找到属于多实验且与该实验相关的图表
    charts = Chart.filter(Chart.project_id == project_id, Chart.name.in_(tag_names))

    db = connect()
    with db.atomic():
        # 必须先清除数据库中的实验数据
        Experiment.delete().where(Experiment.id == experiment_id).execute()
        # 图表无 source 的需要删除
        del_list = []
        for chart in charts:
            if chart.sources.count() == 0:
                del_list.append(chart.id)
        Chart.delete().where(Chart.id.in_(del_list)).execute()
    db.commit()

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
    experiment: Experiment = Experiment.get(experiment_id)
    experiment.update_status(Experiment.STOPPED_STATUS)

    return SUCCESS_200(
        {
            "id": experiment_id,
            "status": Experiment.STOPPED_STATUS,
            "finish_time": create_time(),
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


# 修改实验是否可见
def change_experiment_visibility(experiment_id: int, show: bool):
    """修改实验是否可见

    Parameters
    ----------
    experiment_id : int
        实验id
    show : bool
        在多实验对比图表中，该实验是否可见，true -> 1 , false -> 0

    Returns
    -------
    experiment : dict
        当前实验信息
    """

    try:
        experiment = Experiment.get_by_id(experiment_id)
    except NotExistedError:
        return NOT_FOUND_404("Experiment with id {} does not exist.".format(experiment_id))

    if show:
        experiment.show = 1
    else:
        experiment.show = 0
    experiment.save()
    return SUCCESS_200({"experiment": experiment.__dict__()})
