#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/18 15:28
@File: test_env.py
@IDE: pycharm
@Description:
    测试swanlab.env模块
"""
import os

import pytest

import swanlab
from swanlab.env import SwanLabEnv, is_interactive, get_save_dir


def test_set_default():
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


def test_set_by_netrc():
    """
    测试通过netrc文件设置环境变量
    """
    del os.environ[SwanLabEnv.WEB_HOST.value]
    del os.environ[SwanLabEnv.API_HOST.value]
    del os.environ[SwanLabEnv.RUNTIME.value]
    netrc_path = os.path.join(get_save_dir(), ".netrc")
    with open(netrc_path, "w") as f:
        f.write("machine https://example.cn\nlogin test\npassword 123")
    swanlab.env.SwanLabEnv.set_default()
    assert swanlab.package.get_host_web() == "https://example.cn"
    assert swanlab.package.get_host_api() == "https://example.cn/api"
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


def test_is_interactive():
    # 测试时默认返回true
    assert is_interactive() == True
