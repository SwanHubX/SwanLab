#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/14 15:44
@File: utils.py
@IDE: pycharm
@Description:
    测试控制器中的工具
"""
import pytest
from swanlab.server.controller.utils import *
import random


class TestLTTB:
    @pytest.mark.parametrize("data, sample_length", [
        ([{"data": random.random(), "index": i} for i in range(100)], 100),
        ([{"data": random.random(), "index": i} for i in range(1500)], 1500),
        ([{"data": random.random(), "index": i} for i in range(1502)], 1502),
    ])
    def test_no_sample(self, data, sample_length):
        """
        测试数据长度小于等于采样长度
        """
        result = lttb(data)
        assert len(result) == len(data)

    @pytest.mark.parametrize("data, sample_length", [
        ([{"data": random.random(), "index": i} for i in range(1503)], 1500),
        ([{"data": random.random(), "index": i} for i in range(1600)], 1500),
        ([{"data": random.random(), "index": i} for i in range(5000)], 1500),
    ])
    def test_sample(self, data, sample_length):
        """
        测试数据长度大于采样长度
        """
        result = lttb(data)
        assert len(result) == sample_length
