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
from .utils import is_vaild_ip, get_all_ip


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--debug",
    is_flag=True,
    help="Show more logs when use debug mode",
)
# 控制服务发布的ip地址
@click.option(
    "--host",
    default="0.0.0.0",
    help="When shared, swanlab web will run on localhost",
)
# 控制服务发布的端口，默认5092
@click.option(
    "--port",
    "-p",
    default=5092,
    help="The port of swanlab web, default by 5092",
)
def watch(host, debug, port):
    """Run this command to turn on the swanlab service."""
    # ---------------------------------- 日志等级处理 ----------------------------------

    log_level = "info" if debug else "warning"

    # ---------------------------------- 服务地址处理 ----------------------------------

    # 检查地址是否合法
    if not is_vaild_ip(host):
        raise ValueError(f'ip address "{host}" is not vaild')
    # 拿到当前本机可用的所有ip地址
    ipv4 = get_all_ip()
    if host == "0.0.0.0":
        click.echo(f"swanlab is running...")
        click.echo(f"Available on:")
        for ip in ipv4:
            click.echo(f"  - \033[1mhttp://{ip}:{port}\033[0m")
    elif host in ipv4:
        click.echo(f"swanlab running on \033[1mhttp://{host}:{port}\033[0m")
    else:
        raise ValueError(f'ip address "{host}" is not available, please use one of {ipv4}')
    # ---------------------------------- 启动服务 ----------------------------------
    from ..server import app

    # 使用 uvicorn 启动 FastAPI 应用
    uvicorn.run(app, host=host, port=port, log_level=log_level)


if __name__ == "__main__":
    watch()
