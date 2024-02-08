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
import shutil
from fastapi import Request
from urllib.parse import quote
from ...env import get_swanlog_dir
from ...utils import get_a_lock, COLOR_LIST
from ...utils.file import check_desc_format
import yaml
from typing import List
from ...db.utils.chart import transform_to_multi_exp_charts

from ...db import NotExistedError

# 自定义响应
from ..module.resp import (
    SUCCESS_200,
    DATA_ERROR_500,
    CONFLICT_409,
    NOT_FOUND_404,
)

# 功能性函数
from ..settings import (
    get_tag_dir,
    get_exp_dir,
    get_config_path,
)

# 数据表模型
from ...db import (
    Project,
    Experiment,
    Tag,
    Chart,
    Namespace,
)

# 将查询结果对象转为列表
__to_list = Project.search2list


# ---------------------------------- 通用 ----------------------------------

# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID
# 实验运行状态
RUNNING_STATUS = Experiment.RUNNING_STATUS


def __clear_fields(target: List[dict], fields: List[str]) -> list[dict]:
    """遍历字典列表清除多个字段

    Parameters
    ----------
    target : List[dict]
        需要处理的列表
    fields : List[str]
        需要删除的字段列表

    Returns
    -------
    list[dict]
        处理后的字典列表
    """

    for item in target:
        for field in fields:
            item.pop(field)

    return target


# ---------------------------------- 路由对应的处理函数 ----------------------------------


# 获取项目信息
def get_project_info(project_id: int = DEFAULT_PROJECT_ID) -> dict:
    """
    1. 获取项目信息
    2. 列出当前项目下的所有实验
    3. 获取 swanlog 目录路径

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
        data["logdir"] = get_swanlog_dir()
        experiments = __to_list(project.experiments)
        for experiment in experiments:
            experiment["experiment_id"] = experiment["id"]
        data["experiments"] = experiments
        for experiment in data["experiments"]:
            # 将其中的 project_id 字段去除
            experiment.pop("project_id")
            # 检查配置文件是否存在
            config_path = get_config_path(experiment["run_id"])
            if os.path.exists(config_path):
                # 加载config字段
                with get_a_lock(config_path) as f:
                    experiment["config"] = yaml.load(f, Loader=yaml.FullLoader)
        # 色盘
        data["colors"] = COLOR_LIST
        return SUCCESS_200(data)
    except Exception as e:
        return DATA_ERROR_500(f"Get list experiments failed: {e}")


# 获取项目总结信息
def get_project_summary(project_id: int = DEFAULT_PROJECT_ID) -> dict:
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
            "id": experiment.id,
            "name": experiment.name,
            "run_id": experiment.run_id,
            "tags": [tag["name"] for tag in __to_list(experiment.tags)],
        }
        for experiment in experiments
    ]
    ids = [item["id"] for item in exprs]

    # 根据 id 列表找到所有的 tag，提出不含重复 tag 名的元组
    tags = Tag.filter(Tag.experiment_id.in_(ids))
    # tag_names不用编码，因为前端需要展示
    tag_names = list(set(tag["name"] for tag in __to_list(tags)))

    # 所有总结数据
    data = {}
    # 第一层循环对应实验层，每次探寻一个实验
    for expr in exprs:
        # 第二层为 tag 层，在实验目录下遍历所有 tag，若存在则获取最后一个次提交的值
        experiment_summaries = {}
        for tag in expr["tags"]:
            tag_path = get_tag_dir(expr["run_id"], quote(tag, safe=""))
            if not os.path.exists(tag_path):
                experiment_summaries[tag] = "TypeError"
                continue
            logs = sorted([item for item in os.listdir(tag_path) if item != "_summary.json"])
            with open(os.path.join(tag_path, logs[-1]), mode="r") as f:
                try:
                    tag_data = ujson.load(f)
                    # str 转化的目的是为了防止有些不合规范的数据导致返回体对象化失败
                    experiment_summaries[tag] = str(tag_data["data"][-1]["data"])
                except Exception as e:
                    print(f"[expr: {expr['name']} - {tag}] --- {e}")
                    continue

        data[expr["name"]] = experiment_summaries

    return SUCCESS_200({"tags": tag_names, "summaries": data})


# 修改项目信息
async def update_project_info(request: Request, project_id: int = DEFAULT_PROJECT_ID) -> dict:
    """修改项目信息
    - 项目名
    - 项目描述

    Parameters
    ----------
    project_id : int, optional
        默认为 DEFAULT_PROJECT_ID

    Returns
    -------
    dict
        新项目信息
        - name
        - description
    """

    body = await request.json()

    # 检查格式
    body["description"] = check_desc_format(body["description"], False)

    project = Project.filter(Project.id == project_id).first()
    dict_project = project.__dict__()

    # 检查名字
    if "name" in dict_project and dict_project["name"] == body["name"]:
        pass
    else:
        project.name = body["name"]

    # 检查描述
    if "description" in dict_project and dict_project["description"] == body["description"]:
        pass
    else:
        project.description = body["description"]

    project.save()
    return SUCCESS_200({"updates": body})


# 删除项目
async def delete_project(project_id: int = DEFAULT_PROJECT_ID):
    """删除项目
    1. 删除实验目录
    2. 删除数据库

    Parameters
    ----------
    project_id : int, optional
        默认为第一个项目

    TODO: 原子操作，同时删除项目和实验
    """

    # 检查是否有正在运行的实验
    running_exp = Experiment.filter(Experiment.project_id == project_id, Experiment.status == RUNNING_STATUS).count()
    if running_exp > 0:
        return CONFLICT_409("Can't delete project since there is experiment running")

    project = Project.filter(Project.id == project_id).first()
    run_ids = [experiment["run_id"] for experiment in __to_list(project.experiments)]

    # 删除实验目录
    for run_id in run_ids:
        exp_path = get_exp_dir(run_id)
        if os.path.isdir(exp_path):
            shutil.rmtree(exp_path)

    # 清空数据库
    project.delete_instance()

    return SUCCESS_200({})


async def get_project_charts(project_id: int = DEFAULT_PROJECT_ID) -> dict:
    """获取多实验对比图表数据,并且考虑往期版本兼容性
    1. 如果当前项目的chart字段为0，先生成多实验对比数据，跳转步骤2
    2. 依据规则获取所有实验的图表数据
    """

    # COMPAT 兼容以前没有多实验对比数据的情况
    try:
        transform_to_multi_exp_charts(project_id)
    except IndexError:  # 已经有多实验对比数据
        pass

    # 获取当前项目下所有的多实验对比表
    # 暂时只请求chart.type为default或者line的图表
    allow_types = ["default", "line"]
    multi_charts = Chart.filter(Chart.project_id == project_id, Chart.type.in_(allow_types))
    _sum = multi_charts.count()
    # 获取图表配置
    charts = []
    for _chart in multi_charts:
        # 多实验图表的 source 中不是 tag_name，而是 experiment_name
        sources = []
        # 单箭头是通过外键反向索引的 chart -> source -> tag -> experiment => experiment_name
        for source in _chart.sources:
            sources.append(source.tag_id.experiment_id.name)
        # 当前chart的error字段
        error = {
            source.tag_id.experiment_id.name: Chart.json2dict(source.error) for source in _chart.sources if source.error
        }
        charts.append(
            {
                "id": _chart.id,
                "name": _chart.name,
                "type": _chart.type,
                "reference": _chart.reference,
                "config": _chart.config,
                "description": _chart.description,
                "source": sources,
                "mutli": True,
                "error": error,
            }
        )
    # 获取命名空间配置
    namespaces = Namespace.filter(Namespace.project_id == project_id)
    namespace_list = __clear_fields(__to_list(namespaces), ["project_id", "experiment_id"])
    # 通过 namespace 获取 display，从而获取该 namespace 下的 charts
    for index, namespace in enumerate(namespaces):
        displays = []
        for display in __to_list(namespace.displays):
            # 只获取默认图表
            if display["chart_id"]["type"] in allow_types:
                displays.append(display["chart_id"]["id"])
        if len(displays) > 0:
            namespace_list[index]["charts"] = displays

    # 过滤namespace中的空charts的namespace
    namespace_list = [namespace for namespace in namespace_list if "charts" in namespace]
    # 获取项目下所有实验的图表数据
    return SUCCESS_200({"_sum": _sum, "charts": charts, "namespaces": namespace_list})
