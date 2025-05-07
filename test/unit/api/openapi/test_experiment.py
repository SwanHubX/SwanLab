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


