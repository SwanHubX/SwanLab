#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/16 20:42
@File: create_error_chart.py
@IDE: pycharm
@Description:
    创建错误图表
"""
from tutils import open_dev_mode
import swanlab

swanlab.login(api_key=open_dev_mode())

swanlab.init(
    description="此实验的图表将出现问题",
    log_level="debug",
    mode="cloud",
)

for i in range(10):
    swanlab.log({"success": i, "error": "abc"})

for i in range(10):
    if i % 2 == 1:
        swanlab.log({"nan_test": "nan", "inf_test": "inf"})
    else:
        swanlab.log({"nan_test": i, "inf_test": i})
