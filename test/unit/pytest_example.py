#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-03 00:08:50
@File: test/unit/pytest_example.py
@IDE: vscode
@Description:
    测试用例示例，所有测试用例文件以pytest_开头
"""


def add(x, y):
    """简单的加法函数"""
    return x + y


def test_addition():
    """测试加法函数"""
    assert add(3, 5) == 8
    assert add(-1, 1) == 0
    assert add(2.5, 2.5) == 5


def subtract(x, y):
    """简单的减法函数"""
    return x - y


def test_subtraction():
    """测试减法函数"""
    assert subtract(5, 3) == 2
    assert subtract(1, -1) == 2
    assert subtract(2.5, 1) == 1.5
