#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-19 22:09:09
@File: swanlab\server\controller\project.py
@IDE: vscode
@Description:
    项目相关 api 的处理函数
"""

import os
import ujson
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409
from urllib.parse import unquote
from ..settings import get_logs_dir, get_tag_dir

from ...db import (
    Project,
    Experiment,
    Tag,
)

__to_list = Project.search2list


# ---------------------------------- 通用 ----------------------------------

# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID

# ---------------------------------- 路由对应的处理函数 ----------------------------------


# 列出当前项目下的所有实验
def get_experiments_list(project_id: int = 1) -> dict:
    """
    列出当前项目下的所有实验

    Parameters
    ----------
    project_id : int, optional
        项目id, 默认为DEFAULT_PROJECT_ID

    Returns
    -------
    dict:
        列出当前项目下的所有实验
    """

    try:
        project = Project.filter(Project.id == project_id).first()
        data = project.__dict__()
        data["experiments"] = __to_list(project.experiments)
        return SUCCESS_200(data)
    except Exception as e:
        return DATA_ERROR_500(f"Get list experiments failed: {e}")


# 获取项目总结信息
def get_project_summary(project_id: int = 1) -> dict:
    """
    获取项目下所有实验的总结信息

    Parameters
    ----------
    project_id : int, optional
        项目id, 默认为DEFAULT_PROJECT_ID

    Returns
    -------
    dict:
        项目总结信息
    """

    # 查找所有实验，提出 id 列表
    experiments = Experiment.select().where(Experiment.project_id == project_id)
    exprs = [
        {
            "id": experiment["id"],
            "name": experiment["name"],
            "run_id": experiment["run_id"],
        }
        for experiment in __to_list(experiments)
    ]
    ids = [item["id"] for item in exprs]
    run_ids = [item["run_id"] for item in exprs]

    # 根据 id 列表找到所有的 tag，提出不含重复 tag 名的元组
    tags = Tag.filter(Tag.experiment_id.in_(ids))
    tag_names = list(set(unquote(tag["name"]) for tag in __to_list(tags)))

    # 所有总结数据
    data = {}
    # 第一层循环对应实验层，每次探寻一个实验
    for expr in exprs:
        logs_path: str = get_logs_dir(expr["run_id"])
        # 列出所有 tag 目录
        tag_list: list = os.listdir(logs_path)
        # 第二层为 tag 层，在实验目录下遍历所有 tag，若存在则获取最后一个次提交的值
        experiment_summaries = {}
        for tag in tag_list:
            tag_path = get_tag_dir(expr["run_id"], tag)
            logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
            with open(os.path.join(tag_path, logs[-1]), mode="r") as f:
                try:
                    tag_data = ujson.load(f)
                    experiment_summaries[unquote(tag)] = tag_data["data"][-1]["data"]
                except Exception as e:
                    print(f"[expr: {expr['name']} - {tag}] --- {e}")
                    continue
        data[expr["name"]] = experiment_summaries
    return SUCCESS_200({"tags": tag_names, "summaries": data})
