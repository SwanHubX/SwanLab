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
from ..env import swc
from ..log import swanlog as swl


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
    help="The host of swanlab web, default by 127.0.0.1",
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

    swc.init(swc.getcwd(), "server")
    swl.init(swc.output)
    # ---------------------------------- 服务地址处理 ----------------------------------

    # 检查地址是否合法
    if not is_vaild_ip(host):
        swl.error(f'Ip address "{host}" is not vaild')
        raise ValueError(f'ip address "{host}" is not vaild')
    # 拿到当前本机可用的所有ip地址
    ipv4 = get_all_ip()
    ips = [f"http://{ip}:{port}" for ip in ipv4]
    if host == "0.0.0.0":
        swl.info(f"SwanLab Experiment Dashboard running...")
        swl.info(f"Available on: \n" + "\n".join(ips))
    elif host in ipv4:
        swl.info(f"SwanLab Experiment Dashboard running on http://{host}:{port}")
    else:
        swl.error(f'Ip address "{host}" is not available, please use one of \n' + "\n".join(ips))
        raise ValueError(f'ip address "{host}" is not available, please use one of {ipv4}')
    # ---------------------------------- 启动服务 ----------------------------------
    from ..server import app

    # 使用 uvicorn 启动 FastAPI 应用
    uvicorn.run(app, host=host, port=port, log_level="critical")


if __name__ == "__main__":
    watch()
