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
from fastapi import APIRouter
from ..module.resp import SUCCESS_200, NOT_FOUND_404
from ...env import swc
import os
import ujson
from ...utils import DEFAULT_COLOR

# from ...utils import create_time
from urllib.parse import unquote  # 转码路径参数
from typing import List, Dict
from ...utils import get_a_lock

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
    with get_a_lock(swc.project, "r") as f:
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
    with get_a_lock(swc.project, "r") as f:
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
    path = os.path.join(swc.root, experiment["name"], "logs")
    experiment["tags"] = __list_subdirectories(path)
    experiment["default_color"] = DEFAULT_COLOR
    return SUCCESS_200(experiment)


@router.get("/{experiment_id}/tag/{tag}")
async def get_tag_data(experiment_id: int, tag: str):
    """获取表单数据

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，使用时需要 URIComponent 解码
    """
    tag = unquote(tag)
    # FIXME: 在此处完成num字段的解析
    # num=None: 返回所有数据, num=10: 返回最新的10条数据, num=-1: 返回最后一条数据
    num = None
    # 在experiments列表中查找对应实验的信息
    try:
        experiment_name = __find_experiment(experiment_id)["name"]
    except KeyError as e:
        return NOT_FOUND_404("experiment not found")
    # ---------------------------------- 前置处理 ----------------------------------
    # 获取tag对应的存储目录
    tag_path: str = os.path.join(swc.root, experiment_name, "logs", tag)
    if not os.path.exists(tag_path):
        return NOT_FOUND_404("tag not found")
    # 获取目录下存储的所有数据
    # 降序排列，最新的数据在最前面
    files: list = os.listdir(tag_path)
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

    # ---------------------------------- tag=-1的情况：返回最后一条数据 ----------------------------------

    # FIXME: 如果tag=-1，返回最后一条数据
    if num == -1:
        pass

    # ---------------------------------- tag=其他正数的情况: 返回最新的num条数据 ----------------------------------

    # ---------------------------------- tag=None的情况：返回所有数据 ----------------------------------

    # 阈值，如果数据量大于阈值，只返回阈值条数据
    threshold = 5000
    # 此时count代表总数据量，接下来按量倒叙读取数据
    if count <= threshold:
        # 读取所有数据
        # tag_json是最后一个文件的数据
        # 按顺序读取其他文件的数据
        tag_data_list = __get_immutable_tags(tag_path, files[:-1])
        # 将数据合并
        for data in tag_data_list:
            tag_data.extend(data)
        tag_data.extend(tag_json["data"])
        # 返回数据
        return SUCCESS_200(data={"sum": len(tag_data), "list": tag_data})
    else:
        # TODO 采样读取数据
        raise NotImplementedError("采样读取数据")


@router.get("/{experiment_id}/status")
async def get_experiment_status(experiment_id: int):
    """获取实验状态

    Parameters
    ----------
    experiment_id : int
        实验唯一id，路径传参
    """
    status = __find_experiment(experiment_id)["status"]
    return SUCCESS_200(data={"status": status})


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
    experiment_path: str = os.path.join(swc.root, __find_experiment(experiment_id)["name"], "logs")
    tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
    summaries = []
    for tag in tags:
        tag_path = os.path.join(experiment_path, tag)
        logs = sorted(os.listdir(tag_path))
        with get_a_lock(os.path.join(tag_path, logs[-1]), mode="r") as f:
            data = ujson.load(f)
            summaries.append([tag, data["data"][-1]["data"]])
    return SUCCESS_200(data={"summaries": summaries})


@router.get("/{experiment_id}/log")
async def get_experiment_log(experiment_id: int, page: int):
    """获取收集到的控制台打印

    Parameters
    ----------
    experiment_id : int
        实验唯一ID
    page : int
        分页页码
    """
    # 获取收集到的日志列表
    console_path: str = os.path.join(swc.root, __find_experiment(experiment_id)["name"], "console")
    consoles = [f for f in os.listdir(console_path)]
    total = len(consoles)
    # 如果 page 超出范围
    if not 1 <= page <= total:
        return NOT_FOUND_404("page index out of range")
    # 排序
    consoles = sorted(consoles, key=lambda x: datetime.strptime(x[:-4], "%Y-%m-%d"), reverse=True)
    file_name = consoles[page - 1]
    with get_a_lock(os.path.join(console_path, file_name), mode="r") as f:
        data = f.read()
    return SUCCESS_200(data={"total": total, "logs": data})
