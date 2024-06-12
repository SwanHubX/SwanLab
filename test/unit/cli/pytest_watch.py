#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/11 22:06
@File: pytest_watch.py
@IDE: pycharm
@Description:
    测试cli的watch命令
"""
from swanlab.cli import cli
from click.testing import CliRunner


# noinspection PyTypeChecker
def test_watch_default():
    """
    测试watch命令的默认情况
    """
    runner = CliRunner()
    result = runner.invoke(cli, ["watch"])
    pass
