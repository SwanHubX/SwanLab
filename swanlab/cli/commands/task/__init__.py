#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/17 17:16
@File: __init__.py.py
@IDE: pycharm
@Description:
    启动！
    beta版
"""
from .launch import launch
from .list import list
import click

__all__ = ["task"]


@click.group()
def task():
    pass


# noinspection PyTypeChecker
task.add_command(launch)
# noinspection PyTypeChecker
task.add_command(list)
