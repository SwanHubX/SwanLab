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


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_exp_state():
    """
    获取一个实验的状态
    """
    api = OpenApi()
    res = api.get_exp_state(project="test_project", exp_cuid="test_exp_cuid")
    # 仅 404 情况
    assert isinstance(res, dict)
    assert "state" in res or "code" in res
    if "code" in res:
        assert res["code"] == 404
    elif "state" in res:
        assert res["state"] == "FINISHED" or res["state"] == "RUNNING"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_experiment():
    """
    获取一个实验的详细信息
    """
    api = OpenApi()
    exp_cuid = "ph3oj1b9of9dqj8e38jzl"
    res = api.get_experiment(project="istprvdsbzpwmykkmxekb", exp_cuid=exp_cuid)
    # 仅 404 情况
    assert isinstance(res, dict)
    assert "name" in res or "code" in res
    if "code" in res:
        assert res["code"] == 404
    elif "name" in res:
        assert res["cuid"] == exp_cuid
        assert isinstance(res["name"], str)
        assert isinstance(res["description"], str | None)
        assert isinstance(res["state"], str)
        assert isinstance(res["createdAt"], str)
        assert isinstance(res["finishedAt"], str | None)
        assert isinstance(res["profile"], dict | None)
        if res["profile"] is not None:
            assert isinstance(res["profile"].get("config"), dict | None)
            assert isinstance(res["profile"].get("metadata"), dict | None)
            assert isinstance(res["profile"].get("requirements"), str | None)
            assert isinstance(res["profile"].get("conda"), str | None)
