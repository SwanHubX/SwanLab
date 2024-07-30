#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/30 16:13
@File: stop.py
@IDE: pycharm
@Description:
    停止任务
"""
import click
from .utils import login_init_sid, UseTaskHttp, validate_six_char_string
from swanlab.error import ApiError


@click.command()
@click.argument("cuid", type=str, callback=validate_six_char_string)
def stop(cuid):
    """
    Stop a task by cuid
    """
    login_init_sid()
    with UseTaskHttp() as http:
        try:
            http.patch(f"/task/status", {"cuid": cuid, "status": "STOPPED", "msg": "User stopped by sdk"})
        except ApiError as e:
            if e.resp.status_code == 404:
                raise click.BadParameter("Task not found")
    click.echo("Task stopped successfully, there may be a few minutes of delay online.")
