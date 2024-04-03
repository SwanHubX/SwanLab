#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 16:46
@File: pytest_before_all.py.py
@IDE: pycharm
@Description:
    用于在所有测试用例执行前执行的代码
"""

import pytest


@pytest.fixture(scope="module", autouse=True)
def setup_before_tests():
    # 在运行测试之前执行的前置操作
    print("Setup before running tests")
