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
from tutils.config import nanoid, CONFIG, PACKAGE_PATH
from swanlab.utils.package import get_host_api
import pytest


@pytest.mark.asyncio
async def test_login_success():
    """
    测试登录成功
    """
    login_info = await login_by_key(CONFIG.get('api-key'))
    assert not login_info.is_fail
