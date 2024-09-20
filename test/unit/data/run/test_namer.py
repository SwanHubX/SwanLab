#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/9/20 14:42
@File: test_namer.py
@IDE: pycharm
@Description:
    测试命名器、取色器
"""
from swanlab.data.run import namer


def test_name_no_index():
    name = namer.generate_name()
    assert isinstance(name, str)


def test_name_with_index():
    name = namer.generate_name(4)
    assert isinstance(name, str)
    # 极大数
    name = namer.generate_name(999999999)
    assert isinstance(name, str)


def test_color_no_index():
    colors = namer.generate_colors()
    assert len(colors) == 2
    assert isinstance(colors, tuple)


def test_color_with_index():
    colors = namer.generate_colors(4)
    assert len(colors) == 2
    assert isinstance(colors, tuple)
    # 极大数
    colors = namer.generate_colors(999999999)
    assert len(colors) == 2
    assert isinstance(colors, tuple)
