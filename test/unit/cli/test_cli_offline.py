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


def test_offline_command_creates_settings_file():
    """测试offline命令是否正确创建设置文件"""
    runner = CliRunner()
    
    # 获取swanlog目录路径
    env_key = SwanLabEnv.SWANLOG_FOLDER.value
    logdir = os.environ.get(env_key) or os.path.join(os.getcwd(), "swanlog")
    logdir = os.path.abspath(logdir)
    settings_file = os.path.join(logdir, ".swanlab_settings.json")
    
    # 删除已存在的设置文件
    if os.path.exists(settings_file):
        os.remove(settings_file)

    # 执行offline命令
    result = runner.invoke(cli, ['offline'])

    # 验证命令执行成功
    assert result.exit_code == 0
    assert "✅ SwanLab mode set to offline" in result.output
    assert "📁 Settings file created:" in result.output

    # 验证设置文件被正确创建
    assert os.path.exists(settings_file)
    
    # 验证设置文件内容
    import json
    with open(settings_file, "r", encoding="utf-8") as f:
        settings = json.load(f)
    assert settings.get("mode") == "offline"


def test_offline_command_overwrites_existing_settings():
    """测试offline命令会覆盖已存在的设置文件"""
    runner = CliRunner()
    
    # 获取swanlog目录路径
    env_key = SwanLabEnv.SWANLOG_FOLDER.value
    logdir = os.environ.get(env_key) or os.path.join(os.getcwd(), "swanlog")
    logdir = os.path.abspath(logdir)
    settings_file = os.path.join(logdir, ".swanlab_settings.json")
    
    # 创建一个包含不同模式的设置文件
    import json
    existing_settings = {"mode": "cloud"}
    with open(settings_file, "w", encoding="utf-8") as f:
        json.dump(existing_settings, f)

    # 执行offline命令
    result = runner.invoke(cli, ['offline'])

    # 验证命令执行成功
    assert result.exit_code == 0
    assert "✅ SwanLab mode set to offline" in result.output

    # 验证设置文件被覆盖为offline
    with open(settings_file, "r", encoding="utf-8") as f:
        settings = json.load(f)
    assert settings.get("mode") == "offline"


def test_offline_command_help():
    """测试offline命令的帮助信息"""
    runner = CliRunner()
    result = runner.invoke(cli, ['offline', '--help'])

    # 验证帮助信息包含正确的描述
    assert result.exit_code == 0
    assert "Set SwanLab mode to offline" in result.output
    assert "SWANLAB_MODE" in result.output 