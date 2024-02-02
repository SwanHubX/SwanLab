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


result = json.loads(asyncio.run(get_project_charts()).body)["data"]
assert result["status"] == 200
data = result["data"]


print(result)
