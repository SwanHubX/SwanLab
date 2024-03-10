#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-10 13:53:16
@File: swanlab/server/controller/chart.py
@IDE: vscode
@Description:
    图表相关操作api
"""
from .db import Chart
from .utils import get_exp_charts, get_proj_charts
from ..module import SUCCESS_200, PARAMS_ERROR_422


def update_charts_status(chart_id: int, status: int):
    """
    更新图表状态，status=1表示置顶，status=0表示正常状态，status=-1表示隐藏

    Parameters
    ----------
    chart_id : int
        图表id
    status : int
        状态，1表示置顶，0表示正常状态，-1表示隐藏
    """
    chart: Chart = Chart.get_by_id(chart_id)

    if not chart:
        return PARAMS_ERROR_422("chart_id not exist")
    if status == 1:
        Chart.pin(id=chart_id)
    elif status == -1:
        Chart.hide(id=chart_id)
    else:
        Chart.restore(id=chart_id)
    if chart.project_id:
        chart_list, namespace_list = get_proj_charts(chart.project_id.id)
    else:
        chart_list, namespace_list = get_exp_charts(chart.experiment_id.id)
    for namespace in namespace_list:
        for j, chart_id in enumerate(namespace["charts"]):
            # 在charts中找到对应的chart_id
            namespace["charts"][j] = next((x for x in chart_list if x["id"] == chart_id), None)
    return SUCCESS_200({"groups": namespace_list})
