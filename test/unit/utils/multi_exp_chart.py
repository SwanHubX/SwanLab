#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-01 17:35:21
@File: test/unit/utils/multi_exp_chart.py
@IDE: vscode
@Description:
    测试多实验图表后端数据生成
"""
import os
import sys
import subprocess
from create_exp import run, save_dir
import shutil

if os.path.exists(save_dir):
    shutil.rmtree(save_dir)
os.mkdir(save_dir)


run()

run()
