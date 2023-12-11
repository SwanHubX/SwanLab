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


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--share",
    is_flag=True,
    help="When shared, swanlab web will run on localhost",
)
@click.option(
    "--debug",
    is_flag=True,
    help="Show more logs when use debug mode",
)
@click.option(
    "--port",
    "-p",
    default=5092,
    help="The port of swanlab web, default by 5092",
)
def watch(share, debug, port):
    """Run this command to turn on the swanlab service."""
    # print("share", share)
    # print("debug", debug)
    # print("port", port)
    from ..server import app

    # 服务地址
    host = "localhost" if share else "127.0.0.1"
    # 日志等级
    log_level = "info" if debug else "warning"
    click.echo(f"swanlab running on \033[1mhttp://{host}:{port}\033[0m")

    # 使用 uvicorn 启动 FastAPI 应用
    uvicorn.run(app, host=host, port=port, log_level=log_level)


if __name__ == "__main__":
    cli()
