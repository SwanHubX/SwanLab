"""
@author: cunyue
@file: test_cli_verify.py
@time: 2025/9/23 11:22
@description: 测试cli verify命令
"""

import os

import pytest
from click.testing import CliRunner

import tutils as T
from swanlab.cli.main import cli


# noinspection PyTypeChecker
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_verify_ok():
    runner = CliRunner()
    result = runner.invoke(cli, ["verify"])
    assert result.exit_code == 0
    assert ("You are logged into " + T.WEB_HOST) in result.output


# noinspection PyTypeChecker
def test_verify_error():
    runner = CliRunner()
    if T.SwanLabEnv.API_KEY.value in os.environ:
        del os.environ[T.SwanLabEnv.API_KEY.value]
    result = runner.invoke(cli, ["verify"])
    assert result.exit_code == 1
    assert 'You are not verified. Please use `swanlab login` to login.' in result.output
