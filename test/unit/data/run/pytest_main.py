#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_main.py
@IDE: pycharm
@Description:
    测试SwanLabRun主类
"""
from swanlab.data.run.main import SwanLabRun, get_run


class TestGetRun:

    def test_before_init(self):
        run = get_run()
        assert run is None

    def test_after_init(self):
        run = SwanLabRun()
