#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/13 01:28
@File: pytest_logout.py
@IDE: pycharm
@Description:
    测试登出模块
"""
from click.testing import CliRunner
from swanlab.cli import cli
from tutils import KEY


# noinspection PyTypeChecker
def test_logout_ok(monkeypatch):
    runner = CliRunner()
    # 先登录
    runner.invoke(cli, ["login", "--api-key", KEY])
    monkeypatch.setattr("builtins.input", lambda x: "y")
    result = runner.invoke(cli, ["logout"])
    assert result.exit_code == 0


# noinspection PyTypeChecker
def test_logout_cancel(monkeypatch):
    runner = CliRunner()
    # 先登录
    runner.invoke(cli, ["login", "--api-key", KEY])
    monkeypatch.setattr("builtins.input", lambda x: "n")
    result = runner.invoke(cli, ["logout"])
    assert result.exit_code == 0


# noinspection PyTypeChecker
def test_logout_no_login():
    runner = CliRunner()
    result = runner.invoke(cli, ["logout"])
    assert result.exit_code == 1
    assert "You are not logged in." in result.output
