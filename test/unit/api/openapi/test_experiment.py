#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/7 09:47
@File: test_experiment.py
@IDE: pycharm
@Description:
    测试开放API的实验相关接口
"""

import pytest

import tutils as T
from swanlab import OpenApi
from swanlab.api.openapi.types import ApiResponse, Experiment, Pagination


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_exp_state():
    """
    获取一个实验的状态
    """
    api = OpenApi()
    res = api.get_exp_state(project="test_project", exp_cuid="test_cuid")
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert res.state == "FINISHED" or res.state == "RUNNING"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_experiment():
    """
    获取一个实验的详细信息
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    res = api.get_experiment(project="test_project", exp_cuid=exp_cuid)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, Experiment)

@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_list_project_exps():
    """
    获取一个项目下的实验列表
    """
    api = OpenApi()
    res = api.list_project_exps(project="SwanLab")
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, list)
        for item in res.data:
            assert isinstance(item, Experiment)

@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_exp_summary():
    """
    获取一个实验的Summary信息
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    res = api.get_exp_summary(project="test_project", exp_cuid=exp_cuid)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, dict)
