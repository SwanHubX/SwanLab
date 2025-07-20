#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/1/27 10:00
@File: test_cli_offline.py
@IDE: pycharm
@Description:
    测试offline命令
"""

import os
import pytest
from click.testing import CliRunner

from swanlab.cli.main import cli
from swanlab.env import SwanLabEnv


def test_offline_command_sets_environment_variable():
    """测试offline命令是否正确设置环境变量"""
    runner = CliRunner()
    mode_key = SwanLabEnv.MODE.value
    
    # 确保环境变量不存在或不是offline
    if mode_key in os.environ:
        del os.environ[mode_key]

    # 执行offline命令
    result = runner.invoke(cli, ['offline'])

    # 验证命令执行成功
    assert result.exit_code == 0
    assert "✅ SwanLab mode set to offline" in result.output
    assert "SWANLAB_MODE=offline" in result.output

    # 验证环境变量被正确设置
    assert os.environ.get(mode_key) == "offline"


def test_offline_command_overwrites_existing_mode():
    """测试offline命令会覆盖已存在的模式设置"""
    runner = CliRunner()
    mode_key = SwanLabEnv.MODE.value
    
    # 设置一个不同的模式
    os.environ[mode_key] = "cloud"

    # 执行offline命令
    result = runner.invoke(cli, ['offline'])

    # 验证命令执行成功
    assert result.exit_code == 0
    assert "✅ SwanLab mode set to offline" in result.output

    # 验证环境变量被覆盖为offline
    assert os.environ.get(mode_key) == "offline"


def test_offline_command_help():
    """测试offline命令的帮助信息"""
    runner = CliRunner()
    result = runner.invoke(cli, ['offline', '--help'])

    # 验证帮助信息包含正确的描述
    assert result.exit_code == 0
    assert "Set SwanLab mode to offline" in result.output
    assert "SWANLAB_MODE" in result.output 