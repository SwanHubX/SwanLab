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
    swanlog.disable_log()
    yield
    swanlog.enable_log()
    get_run().finish()


class TestInitMode:
    """
    测试初始化模式时环境变量的设置
    """

    def test_init_disabled(self):
        run = init(mode="disabled", logdir=generate())
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})  # 不会报错
        a = run.settings.run_dir
        assert not os.path.exists(a)

    def test_init_local(self):
        run = init(mode="local")
        assert os.environ[MODE] == "local"
        run.log({"TestInitMode": 1})  # 不会报错

    def test_init_cloud(self):
        run = init(mode="cloud")
        assert os.environ[MODE] == "cloud"
        run.log({"TestInitMode": 1})  # 不会报错
