#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 00:38:23
@File: swanlab\server\api\project.py
@IDE: vscode
@Description:
    项目相关的api，前缀：/project
"""
import os
from fastapi import APIRouter, Request

from ...utils import get_a_lock, create_time
from ..module.resp import SUCCESS_200, DATA_ERROR_500
from ..module import PT
from swanlab.env import swc
import ujson
from urllib.parse import unquote

router = APIRouter()


# 列出当前项目下的所有实验
@router.get("")
async def _():
    """
    获取项目信息，列出当前项目下的所有实验
    """
    try:
        pt = PT()
        return SUCCESS_200(data=pt.get())
    except Exception:
        return DATA_ERROR_500("project data error")


@router.get("/summaries")
async def summaries(experiment_names: str):
    """获取项目的所有实验的总结信息

    Parameters
    ----------
    experiment_names : str
        实验名 —> 通过`,`拼接起来的字符串

    Returns
    -------
    data : dict
        实验总结数据的集合
    """

    column = []
    data = {}
    # 转为列表
    name_list = experiment_names.split(",")
    for name in name_list:
        experiment_path = os.path.join(swc.root, name, "logs")
        tags = [f for f in os.listdir(experiment_path) if os.path.isdir(os.path.join(experiment_path, f))]
        experiment_summaries = {}
        for tag in tags:
            if not unquote(tag) in column:
                column.append(unquote(tag))
            tag_path = os.path.join(experiment_path, tag)
            logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
            with open(os.path.join(tag_path, logs[-1]), mode="r") as f:
                try:
                    tag_data = ujson.load(f)
                except Exception as e:
                    print(f"[expr: {name} - {tag}] --- {e}")
                    continue
                experiment_summaries[unquote(tag)] = tag_data["data"][-1]["data"]
        data[name] = experiment_summaries
    return SUCCESS_200(data={"tags": column, "summaries": data})


@router.patch("/update")
async def update(request: Request):
    body = await request.json()
    file_path = os.path.join(swc.root, "project.json")
    with open(file_path, "r") as f:
        project = ujson.load(f)
    # 检查名字
    if "name" in project and project["name"] == body["name"]:
        pass
    else:
        project.update({"name": body["name"]})
    # 检查描述
    if "description" in project and project["description"] == body["description"]:
        pass
    else:
        project.update({"description": body["description"]})
    # project["update_time"] = create_time()
    # 写入文件
    with get_a_lock(file_path, "w") as f:
        ujson.dump(project, f, indent=4, ensure_ascii=False)
    return SUCCESS_200({"project": project})
