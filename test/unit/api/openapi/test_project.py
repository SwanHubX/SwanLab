#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/8 16:14
@File: test_project.py
@IDE: VSCode
@Description:
    测试开放API的项目相关接口
"""

import pytest

import tutils as T
from swanlab import OpenApi


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_list_projects():
    """
    测试列出一个 workspace 下的所有项目
    """
    api = OpenApi()
    res = api.list_projects(detail=False)
    assert isinstance(res, dict)
    assert "total" in res
    assert "list" in res
    assert isinstance(res["total"], int)
    assert isinstance(res["list"], list)
    assert len(res["list"]) == res["total"]


