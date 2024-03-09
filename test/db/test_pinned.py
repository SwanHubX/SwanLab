#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-09 22:46:04
@File: test/db/test_pinned.py
@IDE: vscode
@Description:
    手动把某个chart置顶
"""
from swanlab.db import connect
from swanlab.db.models.charts import Chart

connect()
chart_id = 1

Chart.pinned(chart_id)
