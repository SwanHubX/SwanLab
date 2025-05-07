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
    # 用户传入 username, workspace, exp_cuid，其中username默认为当前验证登陆的用户，如需看其他团队的实验，可传入对应的username
    res = api.get_exp_state(workspace="test", exp_cuid="test")
    assert isinstance(res, dict)
    assert "state" in res or "code" in res
    if "code" in res:
        assert res["code"] == 404
    elif "state" in res:
        assert res["state"] == "FINISHED" or res["state"] == "RUNNING"


