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
from .status import status
import click

__all__ = ["job"]


@click.group()
def job():
    pass


# noinspection PyTypeChecker
job.add_command(launch)
# noinspection PyTypeChecker
job.add_command(status)
