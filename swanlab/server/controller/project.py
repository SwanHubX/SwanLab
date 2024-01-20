#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-19 22:09:09
@File: swanlab\server\controller\project.py
@IDE: vscode
@Description:
    项目相关 api 的处理函数
"""

from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409
from ...db import (
    Project,
    Experiment,
    Tag,
)

__to_dict = Project.search2dict
__to_list = Project.search2list


# ---------------------------------- 通用 ----------------------------------

# 默认项目 id
DEFAULT_PROJECT_ID = 1
# 默认项目的实例
default_project = Project().filter(Project.id == DEFAULT_PROJECT_ID).first()

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
        data["experiments"] = __to_dict(project.experiments)
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

    # 表头数据
    column = []
    # 总结数据
    data = []

    experiments = Experiment.select().where(Experiment.project_id == project_id)
    experiment_ids = [experiment["id"] for experiment in __to_list(experiments)]
    tags = Tag.filter(Tag.experiment_id.in_(experiment_ids))

    return SUCCESS_200({"tags": __to_dict(tags), "summaries": data})
