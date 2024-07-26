#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/26 17:22
@File: detail.py
@IDE: pycharm
@Description:
    根据cuid获取任务详情
"""
import click
from swanlab.api import get_http
from .utils import TaskModel, login_init_sid
from rich.syntax import Syntax, Console
import json


def validate_six_char_string(_, __, value):
    if value is None:
        raise click.BadParameter('Parameter is required')
    if not isinstance(value, str):
        raise click.BadParameter('Value must be a string')
    if len(value) != 6:
        raise click.BadParameter('String must be exactly 6 characters long')
    return value


@click.command()
@click.argument("cuid", type=str, callback=validate_six_char_string)
def search(cuid):
    """
    Get task detail by cuid
    """
    login_init_sid()
    http = get_http()
    data = http.get(f"/task/{cuid}")
    console = Console()
    s = Syntax(json.dumps(data, indent=4, ensure_ascii=False), "json")
    console.print(s)
