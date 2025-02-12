#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/13 01:11
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试命令登录
"""
import os

import nanoid
import pytest
from click.testing import CliRunner

import tutils as T
from swanlab.cli.main import cli
from swanlab.env import SwanLabEnv
from swanlab.error import ValidationError, APIKeyFormatError
from swanlab.package import get_key


# noinspection PyTypeChecker
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_login_ok():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", T.API_KEY])
    assert result.exit_code == 0
    assert get_key() == T.API_KEY


# noinspection PyTypeChecker
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_login_fail():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", nanoid.generate(alphabet='12345', size=21)])
    assert result.exit_code == 1
    assert isinstance(result.exception, ValidationError)
    result = runner.invoke(cli, ["login", "--api-key", "wrong-key"])
    assert result.exit_code == 1
    assert isinstance(result.exception, APIKeyFormatError)


# noinspection PyTypeChecker
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_login_host():
    """
    测试登录时指定host
    """
    runner = CliRunner()
    del os.environ[SwanLabEnv.API_HOST.value]  # 删除环境变量
    del os.environ[SwanLabEnv.WEB_HOST.value]  # 删除环境变量
    result = runner.invoke(cli, ["login", "--api-key", T.API_KEY, "--host", T.API_HOST.rstrip("/api")])
    assert result.exit_code == 0
    del os.environ[SwanLabEnv.API_HOST.value]  # 删除环境变量
    del os.environ[SwanLabEnv.WEB_HOST.value]  # 删除环境变量
    result = runner.invoke(cli, ["login", "--api-key", T.API_KEY, "--host", "http://wrong-host"])
    assert result.exit_code == 2
