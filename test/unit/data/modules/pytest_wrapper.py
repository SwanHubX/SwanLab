#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/3 15:42
@File: pytest_wrapper.py
@IDE: pycharm
@Description:
    测试包装器
"""
import swanlab.data.modules as M


class TestWrapperLine:
    def test_ok(self):
        """
        正常情况
        """
        data = M.Line(1)
        wrapper = M.DataWrapper("test", [data])
        assert wrapper.is_line is True
