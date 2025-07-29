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
from swanlab.env import SwanLabEnv, is_interactive, get_save_dir, remove_host_suffix


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
        f.write("machine https://example.ai\nlogin test\npassword 123")
    swanlab.env.SwanLabEnv.set_default()
    assert swanlab.package.get_host_web() == "https://example.ai"
    assert swanlab.package.get_host_api() == "https://example.ai/api"
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


# ---------------------------------- 测试移除 host 后缀 ----------------------------------
def test_removes_suffix_when_host_ends_with_suffix():
    host = "example.ai/api"
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == "example.ai"


def test_returns_original_host_when_host_does_not_end_with_suffix():
    host = "example.com"
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == "example.com"


def test_handles_empty_host_and_returns_empty():
    host = ""
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == ""


def test_handles_empty_suffix_and_returns_original_host():
    host = "example.com/api"
    suffix = ""
    assert remove_host_suffix(host, suffix) == "example.com/api"


def test_handles_suffix_longer_than_host_and_returns_original_host():
    host = "api"
    suffix = "example.com/api"
    assert remove_host_suffix(host, suffix) == "api"


def test_handles_host_with_multiple_suffixes():
    host = "example.com/api/v1"
    suffix = "/api/v1"
    assert remove_host_suffix(host, suffix) == "example.com"


def test_handles_host_with_blank_suffix():
    host = "example.com  "
    suffix = "api"
    assert remove_host_suffix(host, suffix) == "example.com"
