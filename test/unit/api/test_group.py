#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/30 21:16
@File: test_group.py
@IDE: pycharm
@Description:
    测试开放API的组织相关接口
"""

import pytest

import tutils as T
from swanlab import OpenApi
from swanlab.api.types import ApiResponse


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_get_workspaces():
    """
    获取用户的所有工作空间
    """
    api = OpenApi()
    r = api.list_workspaces()

    assert isinstance(r, ApiResponse)
    assert isinstance(r.data, list)
