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
from .job import launch_job

__all__ = ["launch"]


def launch():
    launch_job()
