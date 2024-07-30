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
from .utils import TaskModel, login_init_sid, UseTaskHttp
from rich.syntax import Console, Syntax
from swanlab.error import ApiError


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
    login_info = login_init_sid()
    with UseTaskHttp() as http:
        try:
            data = http.get(f"/task/{cuid}")
        except ApiError as e:
            if e.resp.status_code == 404:
                raise click.BadParameter("Task not found")
    tm = TaskModel(login_info.username, data)
    """
    任务名称，python版本，入口文件，任务状态，URL，创建时间，执行时间，结束时间，错误信息
    """
    console = Console()
    print("")
    console.print("[bold]Task Info[/bold]")
    console.print(f"[bold]Task Name:[/bold] [yellow]{tm.name}[/yellow]")
    console.print(f"[bold]Python Version:[/bold] [white]{tm.python}[white]")
    console.print(f"[bold]Entry File:[/bold] [white]{tm.index}[white]")
    icon = '✅'
    if tm.status == 'CRASHED':
        icon = '❌'
    elif tm.status == 'STOPPED':
        icon = '🛑'
    elif tm.status != 'COMPLETED':
        icon = '🏃'
    console.print(f"[bold]Status:[/bold] {icon} {tm.status}")
    console.print(f"[bold]Combo:[/bold] [white]{tm.combo}[/white]")
    tm.url is not None and console.print(f"[bold]SwanLab URL:[/bold] {tm.url}")
    console.print(f"[bold]Created At:[/bold] {tm.created_at}")
    tm.started_at is not None and console.print(f"[bold]Started At:[/bold] {tm.started_at}")
    tm.finished_at is not None and console.print(f"[bold]Finished At:[/bold] {tm.finished_at}")
    if tm.status == 'CRASHED':
        console.print(f"[bold][red]Task Error[/red]:[/bold]\n")
        console.print(Syntax(tm.msg, 'python'))
    print("")  # 加一行空行，与开头一致
