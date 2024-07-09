#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:43
@File: pytest_line.py
@IDE: pycharm
@Description:
    测试折线图模块
"""
import math
import pytest
from swanlab.error import DataTypeError
from swanlab.data.modules import Line


def test_line_ok():
    """
    正常的Line解析
    """
    line = Line(1)
    data, raw = line.parse()
    assert data == 1
    assert raw is None
    assert line.get_chart() == line.Chart.LINE
    assert line.get_section() == "default"
    assert line.get_more() is None
    line = Line("1")
    data, raw = line.parse()
    assert data == 1
    assert raw is None

    class MockType:
        def __float__(self):
            return 1.

    line = Line(MockType())
    data, raw = line.parse()
    assert data == 1
    assert raw is None


def test_line_nan():
    """
    NaN解析
    """
    line = Line("NaN")
    data, raw = line.parse()
    assert data == Line.nan
    assert raw is None
    line = Line(math.nan)
    data, raw = line.parse()
    assert data == Line.nan
    assert raw is None


def test_line_inf():
    """
    INF解析
    """
    line = Line("INF")
    data, raw = line.parse()
    assert data == Line.inf
    assert raw is None
    line = Line(math.inf)
    data, raw = line.parse()
    assert data == Line.inf
    assert raw is None


def test_line_error():
    """
    错误解析
    """
    line = Line("a")
    with pytest.raises(DataTypeError) as e:
        line.parse()
    assert e.value.expected == "float"
    assert e.value.got == "str"
    line = Line([1])
    with pytest.raises(DataTypeError) as e:
        line.parse()
    assert e.value.expected == "float"
    assert e.value.got == "list"
    line = Line(None)
    with pytest.raises(DataTypeError) as e:
        line.parse()
    assert e.value.expected == "float"
    assert e.value.got == "NoneType"
