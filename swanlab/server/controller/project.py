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
from ...db import Project

search_to_dict = Project.search2dict


# 列出当前项目下的所有实验
def list_experiments(project_id: int = 1) -> dict:
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
        return SUCCESS_200({"experiments": search_to_dict(project.experiments)})
    except Exception as e:
        return DATA_ERROR_500(f"Get list experiments failed: {e}")
