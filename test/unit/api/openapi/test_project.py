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
from swanlab.api.openapi.types import ApiResponse, Project, Pagination


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_list_projects():
    """
    测试列出一个 workspace 下的所有项目
    """
    api = OpenApi()
    resp = api.list_projects(detail=False)
    assert isinstance(resp, ApiResponse)
    if resp.code == 200:
        assert isinstance(resp.data, list)
        for item in resp.data:
            assert isinstance(item, Project)
