#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 01:39:04
@File: swanlab\cli\main.py
@IDE: vscode
@Description:
    swanlab脚本命令的主入口
"""

import click
from .utils import is_vaild_ip, URL
from ..utils import FONT
import time
import sys
import atexit


@click.group()
def cli():
    pass


@cli.command()
# 控制服务发布的ip地址
@click.option(
    "--host",
    "-h",
    default="127.0.0.1",
    type=str,
    help="The host of swanlab web, default by 127.0.0.1",
    callback=is_vaild_ip,
)
# 控制服务发布的端口，默认5092
@click.option(
    "--port",
    "-p",
    default=5092,
    type=int,
    help="The port of swanlab web, default by 5092",
)
# 日志等级
@click.option(
    "--log-level",
    default="info",
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    help="The level of log, default by info; You can choose one of [debug, info, warning, error, critical]",
)
def watch(log_level: str, host: tuple, port: int):
    """Run this command to turn on the swanlab service."""
    start = time.time()
    # 导入必要的模块
    from ..log import swanlog as swl
    from ..server import app
    import uvicorn

    # ---------------------------------- 日志等级处理 ----------------------------------
    swl.setLevel(log_level)
    # ---------------------------------- 服务地址处理 ----------------------------------
    # 拿到当前本机可用的所有ip地址
    ip, ipv4 = host
    # ---------------------------------- 日志打印 ----------------------------------
    # 耗时
    take_time = int((time.time() - start) * 1000).__str__() + "ms\n\n"
    # 可用URL
    if URL.is_zero_ip(ip):
        tip = "\n".join([URL(i, port).__str__() for i in ipv4])
    else:
        tip = URL(ip, port).__str__()
    tip = tip + "\n"
    swl.info(f"SwanLab Experiment Dashboard ready in " + FONT.bold(take_time) + tip)

    # ---------------------------------- 启动服务 ----------------------------------
    # 使用 uvicorn 启动 FastAPI 应用，关闭原生日志
    # 使用try expcept捕获退出，涉及端口占用等
    try:
        uvicorn.run(app, host=ip, port=port, log_level="critical")
    except SystemExit as e:
        code = e.code
        if code == 1:
            swl.critical("Error while attempting to bind on address ({}, {}): address already in use".format(ip, port))
        else:
            swl.critical("Unhandled Exit Code: {}".format(code))


if __name__ == "__main__":
    cli()
