#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/2 23:36
@File: test_base.py
@IDE: pycharm
@Description:
    测试基类
"""
from swanlab.data.modules.base import BaseType


class MockType(BaseType):
    def __init__(self, value):
        super().__init__()
        self.value = value

    def parse(self):
        return self.value, None

    def get_chart(self):
        return self.Chart.LINE

    def get_section(self):
        return "MockType"

    def get_more(self):
        return {"value": self.value}

    def get_config(self):
        return {"value": self.value}


def test_base():
    data = MockType(1)
    assert data.parse() == (1, None)
    assert data.get_chart() == data.Chart.LINE
    assert data.get_section() == "MockType"
    assert data.get_more() == {"value": 1}
    assert data.get_config() == {"value": 1}
