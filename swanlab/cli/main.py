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
import uvicorn
from ..server import app

port = 5092


@click.group()
def cli():
    """入口函数"""
    pass


@cli.command()
@click.option(
    "--share",
    default=False,
    help="When shared, swanlab web will run on localhost",
)
@click.option(
    "--debug",
    is_flag=False,
    help="Show more logs when use debug mode",
)
def watch(share, debug):
    """运行此命令开启swanlab服务"""
    host = "localhost" if share else "127.0.0.1"
    log_level = "info" if debug else "warning"
    click.echo(f"swanlab running on \033[1mhttp://{host}:{port}\033[0m")

    # 使用 uvicorn 启动 FastAPI 应用
    uvicorn.run(app, host=host, port=port, log_level=log_level)
