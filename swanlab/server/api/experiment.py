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
from urllib.parse import unquote  # 转码路径参数

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
    experiments: list = ujson.load(open(CONFIG_PATH, "r", encoding="utf-8"))["experiments"]
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
    project: dict = ujson.load(open(CONFIG_PATH, "r", encoding="utf-8"))
    experiments: list = project["experiments"]
    # 在experiments列表中查找对应实验的信息
    for index, experiment in enumerate(experiments):
        if not experiment["experiment_id"] == experiment_id:
            continue
        if experiment["status"] == 0:
            return ResponseBody(500, message="experiment not finished")
        if new_info.get("name") is not None:
            # 如果名字和之前一样
            if new_info["name"] == experiment["name"]:
                return ResponseBody(500, message="experiment name is the same")
            # 检查名字是否重复
            if any(item["name"] == new_info["name"] for item in experiments):
                return ResponseBody(500, message="experiment name already exists")
            # TODO 检查名字是否违规
            # TODO 检查名字的步骤应该封装出来
            # 修改实验目录的名字
            os.rename(
                os.path.join(SWANLAB_LOGS_FOLDER, experiment["name"]),
                os.path.join(SWANLAB_LOGS_FOLDER, new_info["name"]),
            )
            # 在配置中修改
            experiments[index]["name"] = new_info["name"]
        if new_info.get("description") is not None:
            experiment["description"] = new_info["description"]
        project["experiments"] = experiments
        project["update_time"] = create_time()
        ujson.dump(project, open(CONFIG_PATH, "w", encoding="utf-8"), ensure_ascii=False, indent=4)
        return ResponseBody(0)
    return ResponseBody(500, message="experiment not found")


# 获取表单数据
@router.get("/{experiment_id}/{tag}")
async def _(experiment_id: int, tag: str):
    """获取表单数据

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，使用时需要 URIComponent 解码
    """
    # URIComponent 解码 tag
    print(unquote(tag))
    # 读取实验信息内容
    experiments: list = ujson.load(open(CONFIG_PATH, "r", encoding="utf-8"))["experiments"]
    # 在experiments列表中查找对应实验的信息
    tag_info: dict = {}
    for index, experiment in enumerate(experiments):
        if experiment["experiment_id"] == experiment_id:
            # 看看tag是否存在
            indices = [index for index, item in enumerate(experiment["tags"]) if item.get("tag") == unquote(tag)]
            if len(indices) == 1:
                # 存在，保存一下tag的基础信息后退出循环
                tag_info: dict = experiment["tags"][indices[0]]
                tag_info["experiment"]: dict = experiment
                break
            else:
                return ResponseBody(500, message="tag not found")
        # 找到最后一个还不存在，报错
        if index == len(experiments) - 1:
            return ResponseBody(500, message="experiment not found")
    # 获取tag对应的存储目录
    tag_path: str = os.path.join(SWANLAB_LOGS_FOLDER, tag_info.get("experiment").get("name"), tag)
    # 获取目录下存储的所有数据
    files: list = os.listdir(tag_path)
    tag_data: list = []
    for file in files:
        tag_data.extend(ujson.load(open(os.path.join(tag_path, file), "r", encoding="utf-8"))["data"])
    # 返回数据
    return ResponseBody(0, data={"sum": tag_info.get("sum"), "list": tag_data})
