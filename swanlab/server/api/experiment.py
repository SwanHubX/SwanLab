#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 00:39:18
@File: swanlab\server\api\experiment.py
@IDE: vscode
@Description:
    实验相关api，前缀：/experiment
"""
from fastapi import APIRouter
from ..utils import ResponseBody
from ...env import SWANLAB_LOGS_FOLDER
import os
import ujson
from ...utils import create_time

router = APIRouter()
CONFIG_PATH = os.path.join(SWANLAB_LOGS_FOLDER, "project.json")


# 获取当前实验信息
@router.get("/{experiment_id}")
async def _(experiment_id: int):
    """获取当前实验的信息

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    """
    # 读取 project.json 文件内容
    experiments: list = ujson.load(open(CONFIG_PATH, "r"))["experiments"]
    # 在experiments列表中查找对应实验的信息
    for experiment in experiments:
        if experiment["experiment_id"] == experiment_id:
            return ResponseBody(0, data=experiment)
    return ResponseBody(500, message="experiment not found")


# 修改实验的信息：名称/描述
@router.patch("/{experiment_id}")
async def _(experiment_id: int, new_info: dict):
    """修改当前实验的信息

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    new_info: dict
        需要修改的内容
        name: str
        description: str
    """
    # 读取 project.json 文件内容
    project: dict = ujson.load(open(CONFIG_PATH, "r"))
    experiments: list = project["experiments"]
    # 在experiments列表中查找对应实验的信息
    for index, experiment in enumerate(experiments):
        if not experiment["experiment_id"] == experiment_id:
            continue
        if experiment["status"] == 0:
            return ResponseBody(500, message="experiment not finished")
        if "name" in new_info:
            # 如果名字和之前一样
            if new_info["name"] == experiment["name"]:
                return ResponseBody(500, message="experiment name is the same")
            # TODO 检查名字是否重复
            # 修改实验目录的名字
            os.rename(
                os.path.join(SWANLAB_LOGS_FOLDER, experiment["name"]),
                os.path.join(SWANLAB_LOGS_FOLDER, new_info["name"]),
            )
            # 在配置中修改
            experiments[index]["name"] = new_info["name"]
        if "description" in new_info:
            experiment["description"] = new_info["description"]
        project["experiments"] = experiments
        project["update_time"] = create_time()
        ujson.dump(project, open(CONFIG_PATH, "w"))
        return ResponseBody(0, message="experiment updated")
    return ResponseBody(500, message="experiment not found")
