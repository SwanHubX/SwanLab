#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/18 15:28
@File: test_env.py
@IDE: pycharm
@Description:
    测试swanlab.env模块
"""
import pytest

from swanlab.env import SwanLabEnv
import swanlab
import os


def test_default():
    """
    测试获取默认的环境变量
    """
    del os.environ[SwanLabEnv.WEB_HOST.value]
    del os.environ[SwanLabEnv.API_HOST.value]
    del os.environ[SwanLabEnv.RUNTIME.value]
    swanlab.env.SwanLabEnv.set_default()
    assert swanlab.package.get_host_web() == "https://swanlab.cn"
    assert swanlab.package.get_host_api() == "https://api.swanlab.cn/api"
    assert os.getenv(SwanLabEnv.RUNTIME.value) == "user"


def test_check():
    """
    测试检查环境变量
    """
    os.environ[SwanLabEnv.MODE.value] = "124345"
    with pytest.raises(ValueError):
        SwanLabEnv.check()
    os.environ[SwanLabEnv.RUNTIME.value] = "124"
    with pytest.raises(ValueError):
        SwanLabEnv.check()
