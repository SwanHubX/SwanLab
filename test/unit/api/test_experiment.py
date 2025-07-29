#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/7 09:47
@File: test_experiment.py
@IDE: pycharm
@Description:
    测试开放API的实验相关接口
"""

import pandas as pd
import pytest

import tutils as T
from swanlab import OpenApi
from swanlab.api.types import ApiResponse, Experiment


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_experiment():
    """
    获取一个实验的详细信息
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    res = api.get_experiment(project="test_project", exp_id=exp_cuid)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, Experiment)


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_list_experiments():
    """
    获取一个项目下的实验列表
    """
    api = OpenApi()
    res = api.list_experiments(project="SwanLab")
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, list)
        for item in res.data:
            assert isinstance(item, Experiment)


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_summary():
    """
    获取一个实验的Summary信息
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    res = api.get_summary(project="test_project", exp_id=exp_cuid)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, dict)


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_delete_experiment():
    """
    删除一个实验
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    res = api.delete_experiment(project="test_project", exp_id=exp_cuid)
    assert isinstance(res, ApiResponse)
    assert res.code in [204, 404]


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_metrics():
    """
    获取实验的指标数据
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    keys = ["accuracy", "loss"]
    res = api.get_metrics(exp_id=exp_cuid, keys=keys)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, pd.DataFrame)
        assert not res.data.empty
        assert all(key in res.data.columns for key in keys)


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_metrics_duplicate_keys():
    """
    获取实验的指标数据
    """
    api = OpenApi()
    exp_cuid = "test_cuid"
    keys = ["accuracy", "loss", "loss"]
    res = api.get_metrics(exp_id=exp_cuid, keys=keys)
    assert isinstance(res, ApiResponse)
    if res.code == 200:
        assert isinstance(res.data, pd.DataFrame)
        assert not res.data.empty
        assert len(res.data.columns.to_list()) == 2
        assert all(key in res.data.columns for key in ["accuracy", "loss"])
