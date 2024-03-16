#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-09 21:43:59
@File: swanlab/server/controller/utils/charts.py
@IDE: vscode
@Description:
    图表相关函数
"""
from ..db import Chart, Display, Namespace
from typing import List, Union


def get_exp_charts(id: int):
    """
    获取单实验图表数据

    Parameters
    ----------
    id : int
        实验id
    """
    charts: List[Chart] = Chart.filter(Chart.experiment_id == id)
    chart_list = Chart.search2list(charts)
    # 获取每个图表对应的数据源
    for index, chart in enumerate(charts):
        sources = []
        error = {}
        # source->experiment_id
        source_map = {}
        for source in Chart.search2list(chart.sources):
            sources.append(source["tag_id"]["name"])
            source_map[source["tag_id"]["name"]] = id
            if source["error"]:
                error[source["tag_id"]["name"]] = Chart.json2dict(source["error"])
        chart_list[index]["error"] = error
        chart_list[index]["source"] = sources
        chart_list[index]["multi"] = False
        chart_list[index]["source_map"] = source_map

    # 当前实验下的命名空间
    namespaces = Namespace.filter(Namespace.experiment_id == id)
    namespace_list = Namespace.search2list(namespaces)
    # 获取每个命名空间对应的 display
    # display 含有 chart 与 namespace 的对应关系
    for index, namespace in enumerate(namespaces):
        displays = []
        for display in Namespace.search2list(namespace.displays):
            displays.append(display["chart_id"]["id"])
        namespace_list[index]["charts"] = displays
    return get_pinned_and_hidden(chart_list, namespace_list)


def get_proj_charts(id: int):
    """
    获取多实验对比图表数据

    Parameters
    ----------
    id : int
        项目id
    """
    # 获取当前项目下所有的多实验对比表
    # 暂时只请求chart.type为default或者line的图表
    allow_types = ["default", "line", "image", "audio"]
    multi_charts = Chart.filter(Chart.project_id == id, Chart.type.in_(allow_types))
    # 获取图表配置
    charts = []
    # source -> experiment_id
    source_map = {}
    for _chart in multi_charts:
        # 多实验图表的 source 中不是 tag_name，而是 experiment_name
        sources, error = [], {}
        # 单箭头是通过外键反向索引的 chart -> source -> tag -> experiment => experiment_name
        for source in _chart.sources:
            sources.append(source.tag_id.experiment_id.name)
            source_map[source.tag_id.experiment_id.name] = source.tag_id.experiment_id.id
        # 当前chart的error字段，将所有的error字段转换为dict，key为实验名
        for source in _chart.sources:
            if source.error:
                error[source.tag_id.experiment_id.name] = Chart.json2dict(source.error)
        t = _chart.__dict__()
        charts.append({**t, "error": error, "source": sources, "multi": True, "source_map": source_map})
    # 获取命名空间配置
    namespaces = Namespace.filter(Namespace.project_id == id)
    namespace_list = Namespace.search2list(namespaces)
    # 通过 namespace 获取 display，从而获取该 namespace 下的 charts
    for index, namespace in enumerate(namespaces):
        displays = []
        for display in Display.search2list(namespace.displays):
            # 只获取默认图表
            if display["chart_id"]["type"] in allow_types:
                displays.append(display["chart_id"]["id"])
        if len(displays) > 0:
            namespace_list[index]["charts"] = displays
    # 过滤namespace中的空charts的namespace
    namespace_list = [namespace for namespace in namespace_list if "charts" in namespace]

    return get_pinned_and_hidden(charts, namespace_list)


def get_pinned_and_hidden(chart_list: List[dict], namespace_list: List[dict]) -> Union[List[dict], List[dict]]:
    """
    获取置顶图表与隐藏图表
    """
    if not len(chart_list) or not len(namespace_list):
        return chart_list, namespace_list
    # 获取pinned和hidden的开启/关闭状态，通过chart_list的第一个元素的experiment_id或者project_id获取
    first_chart = chart_list[0]
    if first_chart["experiment_id"] is not None:
        exp_or_proj = first_chart["experiment_id"]
    else:
        exp_or_proj = first_chart["project_id"]
    pinned_opened, hidden_opened = exp_or_proj["pinned_opened"], exp_or_proj["hidden_opened"]

    # 遍历chart_list，动态生成pinned与hidden的namespace，这两个namespace的id分别为-1与-2
    pinned_namespace = {
        "id": -1,
        "name": "pinned",
        "charts": [],
        "opened": pinned_opened,
        "experiment_id": first_chart["experiment_id"],
        "project_id": first_chart["project_id"],
    }
    hidden_namespace = {
        "id": -2,
        "name": "hidden",
        "charts": [],
        "opened": hidden_opened,
        "experiment_id": first_chart["experiment_id"],
        "project_id": first_chart["project_id"],
    }
    for chart in chart_list:
        # 如果chart的status为1，则将其加入pinned的namespace，如果是-1加入hidden的namespace
        # 如果是0，则不加入任何namespace
        # 首先将chart对象加入pinned或hidden的namespace，后续会滤除为id
        if chart["status"] == 1:
            pinned_namespace["charts"].append(chart)
        elif chart["status"] == -1:
            hidden_namespace["charts"].append(chart)
        if chart["status"] != 0:
            del_chart_from_namespace(namespace_list, chart["id"])
    # 滤除namespace中的空charts的namespace
    namespace_list = [namespace for namespace in namespace_list if len(namespace["charts"]) != 0]
    # 如果pinned_namespace中有charts，则加入namespace_list的首位
    if len(pinned_namespace["charts"]) > 0:
        # namespaces的charts字段根据每个元素的sort排序，小的在前
        pinned_namespace["charts"].sort(key=lambda x: x["sort"])
        pinned_namespace = {**pinned_namespace, "charts": [chart["id"] for chart in pinned_namespace["charts"]]}
        namespace_list.insert(0, pinned_namespace)
    # 如果hidden_namespace中有charts，则加入namespace_list的末位
    if len(hidden_namespace["charts"]) > 0:
        # namespaces的charts字段根据每个元素的sort排序，小的在前
        hidden_namespace["charts"].sort(key=lambda x: x["sort"])
        hidden_namespace = {**hidden_namespace, "charts": [chart["id"] for chart in hidden_namespace["charts"]]}
        namespace_list.append(hidden_namespace)
    return chart_list, namespace_list


def del_chart_from_namespace(namespace_list: List[dict], chart_id: int) -> List[dict]:
    """
    从命名空间中删除图表
    """
    for namespace in namespace_list:
        if chart_id in namespace["charts"]:
            namespace["charts"].remove(chart_id)
    return namespace_list
