#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/12 21:11
@File: login.py
@IDE: pycharm
@Description:
    登录模块
"""
import os

import click
from swankit.log import FONT

from swanlab.api.auth import terminal_login
from swanlab.env import SwanLabEnv
from swanlab.package import has_api_key, HostFormatter


@click.command()
@click.option(
    "--relogin",
    "-r",
    is_flag=True,
    default=False,
    help="Relogin to the swanlab cloud, it will recover the token file.",
)
@click.option(
    "--api-key",
    "-k",
    default=None,
    type=str,
    help="If you prefer not to engage in commands-line interaction to input the api key, "
    "this will allow automatic login.",
)
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host of the swanlab server.",
)
@click.option(
    "--web-host",
    "-w",
    default=None,
    type=str,
    help="The web host of the swanlab cloud front-end.",
)
def login(api_key: str, relogin: bool, host: str, web_host: str):
    """Login to the swanlab cloud."""
    if not relogin and has_api_key():
        # 此时代表token已经获取，需要打印一条信息：已经登录
        command = FONT.bold("swanlab login --relogin")
        tip = FONT.swanlab("You are already logged in. Use `" + command + "` to force relogin.")
        return print(tip)
    # 清除环境变量
    if relogin:
        del os.environ[SwanLabEnv.API_HOST.value]
        del os.environ[SwanLabEnv.WEB_HOST.value]
    try:
        HostFormatter(host, web_host)()
    except ValueError as e:
        raise click.BadParameter(str(e))
    # 进行登录，此时将直接覆盖本地token文件
    login_info = terminal_login(api_key)
    print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")
