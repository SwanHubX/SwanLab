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
def test_get_workspaces():
    """
    获取用户的所有工作空间
    """
    api = OpenApi()
    # 用户传入 username, workspace, exp_cuid
    res = api.get_exp_state(username="", workspace="", exp_cuid="")
    assert isinstance(res, dict)
    assert "state" in res
