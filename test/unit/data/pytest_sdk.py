#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_sdk.py
@IDE: pycharm
@Description:
    测试sdk的一些api
"""
from swanlab import init
from swanlab.env import (
    reset_env,
    STRICT_MODE,
    MODE,
    ROOT
)
from swanlab.log import swanlog
from swanlab.data.run import get_run
from nanoid import generate
import pytest
import os


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    reset_env()
    root = os.environ[ROOT]
    if STRICT_MODE in os.environ:
        del os.environ[STRICT_MODE]
    if MODE in os.environ:
        del os.environ[MODE]
    del os.environ[ROOT]
    swanlog.disable_log()
    yield
    swanlog.enable_log()
    reset_env()
    if STRICT_MODE in os.environ:
        del os.environ[STRICT_MODE]
    if MODE in os.environ:
        del os.environ[MODE]
    os.environ[ROOT] = root
    get_run().finish()


class TestInitMode:
    """
    测试初始化模式时环境变量的设置
    """

    def test_init_disabled(self):
        run = init(mode="disabled", logdir=generate())
        assert os.environ[MODE] == "disabled"
        run.log({"test": 1})

    def test_init_local(self):
        init(mode="local")
        assert os.environ[MODE] == "local"

    def test_init_cloud(self):
        init(mode="cloud")
        assert os.environ[MODE] == "cloud"
