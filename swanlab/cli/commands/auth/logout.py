#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/12 21:11
@File: logout.py
@IDE: pycharm
@Description:
    登出模块
"""
import shutil
import sys

import click
from rich.text import Text

from swanlab.env import get_save_dir
from swanlab.log import swanlog
from swanlab.package import has_api_key


@click.command()
def logout():
    """Logout to the swanlab cloud."""
    command = Text("swanlab login", "bold")
    if has_api_key():
        # 如果已经是登录状态，那么则询问用户是否确认，如果确认则删除token文件夹

        swanlog.info("Are you sure you want to logout? (y/N): ", end='')
        confirm = input()
        if confirm.lower() == "y":
            try:
                shutil.rmtree(get_save_dir())
                return swanlog.info(" Logout successfully. You can use `", command, "` to login again.", sep="")
            except Exception as e:
                return swanlog.info("Logout failed. Reason:" + str(e))
        else:
            return swanlog.info("Logout canceled.")
    # 如果还未登录，则不做任何处理，并告知用户如何登录
    swanlog.info("You are not logged in. If you want to login in, please use `", command, "` to login.", sep="")
    return sys.exit(1)
