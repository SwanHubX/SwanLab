#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:16
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试登录
"""
from swanlab.api.auth.login import login_by_key
from tutils.config import CONFIG
import pytest


@pytest.mark.asyncio
async def test_login_success():
    """
    测试登录成功
    """
    login_info = await login_by_key(CONFIG.get('api-key'), save=False)
    assert not login_info.is_fail
    assert login_info.api_key == CONFIG.get('api-key')
    assert login_info.__str__() == "Login success"


@pytest.mark.asyncio
async def test_login_error_key():
    """
    测试登录失败, 错误的key
    """
    login_info = await login_by_key("wrong-key", save=False)
    assert login_info.is_fail
    assert login_info.api_key is None
    assert login_info.__str__() == "Error api key"


@pytest.mark.asyncio
async def test_login_no_key():
    """
    测试登录失败, 未输入key
    """
    login_info = await login_by_key(None, save=False)
    assert login_info.is_fail
    assert login_info.api_key is None
    assert login_info.__str__() == "Error api key"
