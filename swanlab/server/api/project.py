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
import shutil
from fastapi import APIRouter, Request

from ...utils import get_a_lock
from ...utils.file import check_desc_format
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409
from swanlab.env import get_swanlog_dir
import ujson
from urllib.parse import unquote
from ..settings import PROJECT_PATH, SWANLOG_DIR

router = APIRouter()


# 列出当前项目下的所有实验
@router.get("")
async def _():
    """
    获取项目信息，列出当前项目下的所有实验
    """
    try:
        return SUCCESS_200({})
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
        experiment_path = os.path.join(get_swanlog_dir(), name, "logs")
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
    """更新项目配置信息

    Parameters
    ----------
    request : Request
        fastapi中请求对象，可以通过上面的 json 方法将请求体解析成 dict 格式

    Returns
    -------
    project：Obj
        返回 project.json 的所有内容，目的是方便前端在修改信息后重置 pinia 的状态
    """
    body = await request.json()
    # 检查格式
    body["description"] = check_desc_format(body["description"], False)

    with open(PROJECT_PATH, "r") as f:
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
    with get_a_lock(PROJECT_PATH, "w") as f:
        ujson.dump(project, f, indent=4, ensure_ascii=False)
    return SUCCESS_200({"project": project})


@router.delete("")
async def delete():
    """删除项目
    1. 如果项目中还有实验正在运行，应不予删除，并报错提示
    2. 仅删除 project.json 中记录的内容，即 swanlab 自己生产的东西
    """

    folders = []
    # 检测是否有正在运行的实验
    with open(PROJECT_PATH, "r") as f:
        project = ujson.load(f)
    for item in project["experiments"]:
        if item["status"] == 0:
            return CONFLICT_409("Can't delete project since there is experiment running")
        else:
            folders.append(item["name"])

    # 清空除了日志以外的项目文件
    # 需要清除的有：1. 实验目录 2. project.json
    for folder in folders:
        folder_path = os.path.join(SWANLOG_DIR, folder)
        if os.path.isdir(folder_path):
            shutil.rmtree(folder_path)
    os.remove(PROJECT_PATH)

    return SUCCESS_200({})
