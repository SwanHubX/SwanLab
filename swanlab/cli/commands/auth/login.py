#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/12 21:11
@File: login.py
@IDE: pycharm
@Description:
    登录模块
"""
import click
from swanlab.package import is_login
from swanlab.api.auth import terminal_login
from swankit.log import FONT


@click.command()
@click.option(
    "--relogin", "-r",
    is_flag=True,
    default=False,
    help="Relogin to the swanlab cloud, it will recover the token file.",
)
@click.option(
    "--api-key", "-k",
    default=None,
    type=str,
    help="If you prefer not to engage in commands-line interaction to input the api key, "
         "this will allow automatic login.",
)
def login(api_key: str, relogin: bool):
    """Login to the swanlab cloud."""
    if not relogin and is_login():
        # 此时代表token已经获取，需要打印一条信息：已经登录
        command = FONT.bold("swanlab login --relogin")
        tip = FONT.swanlab("You are already logged in. Use `" + command + "` to force relogin.")
        return print(tip)
    # 进行登录，此时将直接覆盖本地token文件
    login_info = terminal_login(api_key)
    print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")
