#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/13 01:11
@File: pytest_login.py
@IDE: pycharm
@Description:
    测试命令登录
"""

from swanlab.env import get_swanlab_folder
from swanlab.package import get_key, get_host_api
from click.testing import CliRunner
from swanlab.cli.main import cli
from tutils import KEY
from swanlab.error import ValidationError
import os


# noinspection PyTypeChecker
def test_login_ok():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", KEY])
    assert result.exit_code == 0
    path = os.path.join(get_swanlab_folder(), ".netrc")
    assert os.path.exists(path)
    assert get_key(path, get_host_api())[2] == KEY


# noinspection PyTypeChecker
def test_login_fail():
    runner = CliRunner()
    result = runner.invoke(cli, ["login", "--api-key", "123"])
    assert result.exit_code == 1
    assert isinstance(result.exception, ValidationError)
