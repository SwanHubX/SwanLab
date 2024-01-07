#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 00:39:18
@File: swanlab\server\api\experiment.py
@IDE: vscode
@Description:
    实验相关api，前缀：/experiment
"""
from datetime import datetime
import shutil
from fastapi import APIRouter, Request

from ...utils.file import check_exp_name_format, check_desc_format
from ..module.resp import SUCCESS_200, NOT_FOUND_404, PARAMS_ERROR_422, Conflict_409
import os
import ujson
from urllib.parse import quote, unquote  # 转码路径参数
from typing import List, Dict
from ..settings import PROJECT_PATH, SWANLOG_DIR
from ...utils import get_a_lock, create_time, DEFAULT_COLOR
from ...log import swanlog as swl

router = APIRouter()


# ---------------------------------- 工具函数 ----------------------------------
def __find_experiment(experiment_id: int) -> dict:
    """在实验列表中查找对应id的实验

    Parameters
    ----------
    experiment_id : int
        实验id

    Returns
    -------
    dict
        实验信息
    """
    with get_a_lock(PROJECT_PATH, "r") as f:
        experiments: list = ujson.load(f)["experiments"]
    for experiment in experiments:
        if experiment["experiment_id"] == experiment_id:
            return experiment
    # 还不存在就报错
    raise KeyError(f'experiment id "{experiment_id}" not found')


def __get_immutable_tags(base_path: str, paths: list) -> List[List[Dict]]:
    """批量读取tag数据，你需要保证paths中的tag文件不再更新

    Parameters
    ----------
    base_path : str
        实验tag的存储路径
    paths : list
        tag的存储路径列表
    """
    tag_data_list: List[List[Dict]] = []
    for path in paths:
        # 读取tag数据，由于目前在设计上这些文件不会再被修改，所以不需要加锁
        with open(os.path.join(base_path, path), "r") as f:
            tag_data_list.append(ujson.load(f)["data"])
    return tag_data_list


def __list_subdirectories(folder_path: str) -> List[str]:
    """列出文件夹下的所有子文件夹，过滤子文件

    Parameters
    ----------
    folder_path : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    items = os.listdir(folder_path)

    # 使用列表推导式筛选出所有的子文件夹
    subdirectories = [item for item in items if os.path.isdir(os.path.join(folder_path, item))]

    return subdirectories


def __get_charts(chart_path: str):
    """获取实验图表配置，并且向下兼容，如果配置不完整，自动补全

    Parameters
    ----------
    chart_path : str
        实验图表配置文件路径
    """
    with get_a_lock(chart_path, "r+") as f:
        chart: dict = ujson.load(f)
        # COMPAT 如果chart不存在namespaces且charts有东西，生成它
        compat = not chart.get("namespaces") and len(chart["charts"])
        if compat:
            # 提示用户，配置将更新
            swl.warning(
                "The configuration of the chart is somewhat outdated. SwanLab will automatically make some updates to this configuration."
            )
            # 遍历chart[charts],写入chart_id
            charts = [c["chart_id"] for c in chart["charts"]]
            ns = {"namespace": "default", "charts": charts}
            chart["namespaces"] = [ns]
            # 写入文件
            f.truncate(0)
            f.seek(0)
            ujson.dump(chart, f, ensure_ascii=False)
    return chart


# ---------------------------------- 业务路由 ----------------------------------


@router.get("/{experiment_id}")
async def get_experiment(experiment_id: int):
    """获取当前实验的信息

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    """
    # 读取 project.json 文件内容
    with get_a_lock(PROJECT_PATH, "r") as f:
        experiments: list = ujson.load(f)["experiments"]
    # 在experiments列表中查找对应实验的信息
    experiment = None
    for ex in experiments:
        if ex["experiment_id"] == experiment_id:
            experiment = ex
            break
    # 如果没有找到，即实验不存在
    if experiment is None:
        return NOT_FOUND_404()
    # 生成实验存储路径
    path = os.path.join(SWANLOG_DIR, experiment["name"], "logs")
    experiment["tags"] = __list_subdirectories(path)
    experiment["default_color"] = DEFAULT_COLOR
    return SUCCESS_200(experiment)


# COMPAT 由于fastapi不支持%2F的路径转换，所以采用通配符:path，并且在下面将path进行quote
@router.get("/{experiment_id}/tag/{tag:path}")
async def get_tag_data(experiment_id: int, tag: str):
    """获取表单数据

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，使用时需要 URIComponent 解码
    """
    tag = quote(tag, safe="")
    # 在experiments列表中查找对应实验的信息
    try:
        experiment_name = __find_experiment(experiment_id)["name"]
    except KeyError as e:
        return NOT_FOUND_404("experiment not found")
    # ---------------------------------- 前置处理 ----------------------------------
    # 获取tag对应的存储目录
    tag_path: str = os.path.join(SWANLOG_DIR, experiment_name, "logs", tag)
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
    tag_data_list = __get_immutable_tags(tag_path, files[:-1])
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
        swl.warning(warn)
        # 遍历tag_data，获取最大值和最小值
        # 提取 data 字段的值
        data_values = [entry["data"] for entry in tag_data]
        # 获取最大值和最小值
        max_value = max(data_values)
        min_value = min(data_values)
    return SUCCESS_200(data={"sum": len(tag_data), "max": max_value, "min": min_value, "list": tag_data})


@router.get("/{experiment_id}/status")
async def get_experiment_status(experiment_id: int):
    """获取实验状态以及实验图表配置，用于实时更新实验状态

    Parameters
    ----------
    experiment_id : int
        实验唯一id，路径传参
    """
    exp = __find_experiment(experiment_id)
    status = exp["status"]
    chart_path: str = os.path.join(SWANLOG_DIR, exp["name"], "chart.json")
    charts = __get_charts(chart_path)
    return SUCCESS_200(data={"status": status, "charts": charts})


@router.get("/{experiment_id}/summary")
async def get_experiment_summary(experiment_id: int):
    """获取实验的总结数据——每个tag的最后一个setp的data

    Parameters
    ----------
    experiment_id : int
        实验id

    Returns
    -------
    array
        每个tag的最后一个数据
    """
    experiment_path: str = os.path.join(SWANLOG_DIR, __find_experiment(experiment_id)["name"], "logs")
    tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
    tags = [item for item in tags if item != "_summary.json"]
    summaries = []
    for tag in tags:
        tag_path = os.path.join(experiment_path, tag)
        logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
        with get_a_lock(os.path.join(tag_path, logs[-1]), mode="r") as f:
            data = ujson.load(f)
            # 保留4位有效数字
            data = round(data["data"][-1]["data"], 4)
            summaries.append({"key": unquote(tag), "value": data})
    return SUCCESS_200(data={"summaries": summaries})


MAX_NUM = 6000


@router.get("/{experiment_id}/recent_log")
async def get_recent_experiment_log(experiment_id: int):
    """一下返回最多 MAX_NUM 条打印记录

    Parameters
    ----------
    experiment_id : int
        实验唯一ID
    MAX_NUM : int
        最多返回这么多条
    """
    console_path: str = os.path.join(SWANLOG_DIR, __find_experiment(experiment_id)["name"], "console")
    consoles: list = [f for f in os.listdir(console_path)]
    # 含有error.log，在返回值中带上其中的错误信息
    error = None
    if "error.log" in consoles:
        with open(os.path.join(console_path, "error.log"), mode="r") as f:
            error = f.read().split("\n")
        # 在consoles里删除error.log
        consoles.remove("error.log")
    total: int = len(consoles)
    # 如果 total 大于 1, 按照时间排序
    if total > 1:
        consoles = sorted(consoles, key=lambda x: datetime.strptime(x[:-4], "%Y-%m-%d"), reverse=True)
    logs = []
    # current_page = total
    for index, f in enumerate(consoles, start=1):
        with open(os.path.join(console_path, f), mode="r") as f:
            logs.extend(f.read().split("\n"))
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


@router.get("/{experiment_id}/chart")
async def get_experimet_charts(experiment_id: int):
    chart_path: str = os.path.join(SWANLOG_DIR, __find_experiment(experiment_id)["name"], "chart.json")
    chart = __get_charts(chart_path)
    return SUCCESS_200(chart)


@router.get("/{experiment_id}/stop")
async def stop_experiment(experiment_id: int):
    """停止实验

    Parameters
    ----------
    experiment_id : int
        实验唯一ID
    """
    with open(PROJECT_PATH, mode="r", encoding="utf-8") as f:
        config = ujson.load(f)
    # 获取需要停止的实验在配置中的索引
    index = next((index for index, d in enumerate(config["experiments"]) if d["experiment_id"] == experiment_id), None)
    # 修改对应实验的状态
    if not config["experiments"][index]["status"] == 0:
        # 不在运行中的状态不予修改
        return Exception("Experiment status is not running")
    config["experiments"][index]["status"] = -1
    config["experiments"][index]["update_time"] = create_time()
    with get_a_lock(PROJECT_PATH, "w") as f:
        ujson.dump(config, f, ensure_ascii=False, indent=4)
    return SUCCESS_200({"update_time": create_time()})


@router.patch("/{experiment_id}")
async def update_experiment_config(experiment_id: int, request: Request):
    """修改实验的元信息

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
    object
    """
    body: dict = await request.json()
    # 校验参数
    check_exp_name_format(body["name"], False)
    body["description"] = check_desc_format(body["description"], False)

    with open(PROJECT_PATH, "r") as f:
        project = ujson.load(f)
    experiment = __find_experiment(experiment_id)
    # 寻找实验在列表中对应的 index
    experiment_index = project["experiments"].index(experiment)

    # 修改实验名称
    if not experiment["name"] == body["name"]:
        # 修改实验名称
        if not experiment["name"] == body["name"]:
            # 检测实验名是否重复
            for expr in project["experiments"]:
                if expr["name"] == body["name"]:
                    return PARAMS_ERROR_422("Experiment's target name already exists")
        project["experiments"][experiment_index]["name"] = body["name"]
        # 修改实验目录名
        old_path = os.path.join(SWANLOG_DIR, experiment["name"])
        new_path = os.path.join(SWANLOG_DIR, body["name"])
        os.rename(old_path, new_path)

    # 修改实验描述
    if not experiment["description"] == body["description"]:
        project["experiments"][experiment_index]["description"] = body["description"]
    with get_a_lock(PROJECT_PATH, "w") as f:
        ujson.dump(project, f, indent=4, ensure_ascii=False)

    return SUCCESS_200({"experiment": project["experiments"][experiment_index]})


@router.delete("/{experiment_id}")
async def delete_experiment(experiment_id: int):
    """删除实验

    注意，需要先判断当前实验是否正在运行中，不可删除运行中的实验

    Parameters
    ----------
    experiment_id : Int
        实验唯一ID

    Returns
    -------
    project : Dictionary
        删除实验后的项目信息，提供给前端更新界面
    """
    experiment: dict = __find_experiment(experiment_id)

    # 判断状态，运行中则不可删除
    if experiment["status"] == 0:
        return Conflict_409("Can't delete experiment since experiment is running")

    # 可以删除
    # 1. 删除 project.json 中的实验记录
    with get_a_lock(PROJECT_PATH) as f:
        project: dict
        with open(PROJECT_PATH, "r") as file_read:
            project = ujson.load(file_read)
            # project["_sum"] = project["_sum"] - 1
            project["experiments"].remove(experiment)

        with open(PROJECT_PATH, "w") as file_write:
            ujson.dump(project, file_write, indent=4, ensure_ascii=False)
    # 2. 删除实验目录
    shutil.rmtree(os.path.join(SWANLOG_DIR, experiment["name"]))

    return SUCCESS_200({"project": project})
