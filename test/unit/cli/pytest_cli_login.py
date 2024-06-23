#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/13 01:11
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试命令登录
"""
from swanlab.package import get_key
from click.testing import CliRunner
from swanlab.cli.main import cli
from tutils import TEST_CLOUD_KEY
from swanlab.error import ValidationError
import tutils as T
import pytest


# noinspection PyTypeChecker
@pytest.mark.skipif(T.TEST_CLOUD_SKIP, reason="skip cloud test")
def test_login_ok():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", TEST_CLOUD_KEY])
    assert result.exit_code == 0
    assert get_key() == TEST_CLOUD_KEY


# noinspection PyTypeChecker
@pytest.mark.skipif(T.TEST_CLOUD_SKIP, reason="skip cloud test")
def test_login_fail():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", "123"])
    assert result.exit_code == 1
    assert isinstance(result.exception, ValidationError)
