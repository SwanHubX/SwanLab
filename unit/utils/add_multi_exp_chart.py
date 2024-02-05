#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-03 17:01:09
@File: unit\utils\add_multi_exp_chart.py
@IDE: vscode
@Description:
    添加多实验对比图表
"""

import create_exp
from swanlab.db import *
from swanlab.db.utils.chart import add_multi_chart


# 连接数据库
connect()

res = add_multi_chart(
    project_id=1,
    tag_id=3,
    chart_id=9,
)

print(res)
