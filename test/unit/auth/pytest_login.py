#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:16
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试登录
"""
from swanlab.auth.login import login_by_key
import pytest


@pytest.mark.asyncio
async def test_login_success():
    """
    测试登录成功
    """
    assert 1 == 1
