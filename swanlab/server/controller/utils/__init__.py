#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-09 21:41:37
@File: swanlab/server/controller/utils/__init__.py
@IDE: vscode
@Description:
    共享工具函数
"""
from .charts import get_exp_charts, get_proj_charts
from .tag import read_tag_data, get_tag_files, LOGS_CONFIGS, lttb
from typing import List


def clear_field(target: List[dict], field: str) -> List[dict]:
    """遍历字典列表清除某个字段

    Parameters
    ----------
    target : List[dict]
        需要处理的列表
    field : str
        需要删除的字段

    Returns
    -------
    List[dict]
        处理后的字典列表
    """

    for item in target:
        item.pop(field)

    return target
