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
from .list import list
from .search import search
from .stop import stop
import click

__all__ = ["task"]


@click.group()
def task():
    """
    Beta Function: List, modify, query task information.
    """
    pass


# noinspection PyTypeChecker
task.add_command(list)
# noinspection PyTypeChecker
task.add_command(search)
# noinspection PyTypeChecker
task.add_command(stop)
