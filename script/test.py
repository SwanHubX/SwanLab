#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 16:41
@File: test.py.py
@IDE: pycharm
@Description:
    运行单元测试，运行时应该在当前项目根目录下运行
"""
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.system("pytest test/unit")
