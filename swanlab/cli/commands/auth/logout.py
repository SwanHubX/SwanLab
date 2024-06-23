#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/12 21:11
@File: logout.py
@IDE: pycharm
@Description:
    登出模块
"""
import click
from swanlab.package import is_login
from swankit.log import FONT
from swanlab.env import get_save_dir
import shutil
import sys


@click.command()
def logout():
    """Logout to the swanlab cloud."""
    command = FONT.bold("swanlab login")
    if is_login():
        # 如果已经是登录状态，那么则询问用户是否确认，如果确认则删除token文件夹
        confirm = input(FONT.swanlab("Are you sure you want to logout? (y/N): "))
        if confirm.lower() == "y":
            try:
                shutil.rmtree(get_save_dir())
                return print(FONT.swanlab("Logout successfully. You can use `" + command + "` to login again."))
            except Exception as e:
                return print(FONT.swanlab("Logout failed. Reason:" + str(e)))
        else:
            return print(FONT.swanlab("Logout canceled."))
    # 如果还未登录，则不做任何处理，并告知用户如何登录
    tip = FONT.swanlab("You are not logged in. If you want to login in, please use `" + command + "` to login.")
    print(tip)
    return sys.exit(1)
