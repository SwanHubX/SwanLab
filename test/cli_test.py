#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-17 15:26:00
@File: test/cil_test.py
@IDE: vscode
@Description:
    测试cil模块
"""
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from swanlab.cli.main import watch


if __name__ == "__main__":
    watch()
