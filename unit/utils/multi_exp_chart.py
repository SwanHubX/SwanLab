#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-01 17:35:21
@File: test/unit/utils/multi_exp_chart.py
@IDE: vscode
@Description:
    测试多实验图表后端数据生成
    运行此测试文件之前，请先运行create_exp.py生成实验数据，至少生成两个实验以上的数据
"""
import create_exp
import asyncio
from swanlab.db import *
import json

# 连接数据库
connect()
from swanlab.server.controller.project import get_project_charts

"""
在开发时，先运行接口，完成接口部分（“向下兼容”）的开发，然后完成数据库部分的开发
"""


result = json.loads(asyncio.run(get_project_charts()).body)

# ---------------------------------- 请求状态 ----------------------------------

if result["code"] == 3404:
    print(result["message"])

assert result["code"] == 0

# ---------------------------------- 响应体字段 ----------------------------------

data = result["data"]

assert "_sum" in data, "data中缺少 _sum 字段"
assert "charts" in data, "data中缺少 charts 字段"
assert "namespaces" in data, "data中缺少 namespaces 字段"

# ---------------------------------- 字段类型检测 ----------------------------------

# _sum
assert isinstance(data["_sum"], (int, float)), "_sum 字段的值不是数字"
# charts
charts_value = data["charts"]
assert isinstance(charts_value, list), "charts 字段的值不是列表"
for chart in charts_value:
    assert isinstance(chart, dict), "charts 列表中的元素不是字典"
    assert all(key in chart for key in ("id", "name", "source")), "charts 字典缺少 id, name 或 source 字段"
# namesapces
namespaces_value = data["namespaces"]
assert isinstance(namespaces_value, list), "namespaces 字段的值不是列表"
for namespace in namespaces_value:
    assert isinstance(namespace, dict), "namespaces 列表中的元素不是字典"
    assert all(key in namespace for key in ("id", "name", "charts")), "namespaces 字典缺少 id, name 或 charts 字段"

# ---------------------------------- 字段值检测 ----------------------------------
# 查询当前数据库中的项目图表
charts = Chart.search2list(Chart.filter(project_id=Project.DEFAULT_PROJECT_ID))
# 这些图表应该存在于data.charts中
for chart in charts:
    assert any(c["id"] == chart["id"] for c in charts_value), f"图表 {chart.id} 不在 data.charts 中"
# 查询当前数据库中的项目命名空间
namespaces = Namespace.search2list(Namespace.filter(project_id=Project.DEFAULT_PROJECT_ID))
# 这些命名空间应该存在于data.namespaces中
for namespace in namespaces:
    assert any(n["id"] == namespace["id"] for n in namespaces_value), f"命名空间 {namespace.id} 不在 data.namespaces 中"
# 这些命名空间下的图表应该存在于data.namespaces.charts中，并且与数据库中的一致
for namespace in namespaces_value:
    namespace_id = namespace["id"]
    namespace_charts = Display.search2list(Display.filter(namespace_id=namespace_id))
    for chart in namespace_charts:
        assert any(
            c == chart["chart_id"]["id"] for c in namespace["charts"]
        ), f"图表 {chart.id} 不在 data.namespaces.charts 中"
print(data)
