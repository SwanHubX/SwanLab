#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 15:30:24
@File: swanlab/auth/login.py
@IDE: vscode
@Description:
    用户登录
"""

import asyncio
import requests
from ..utils import FONT, get_package_version
import sys


async def _login(api_key: str):
    """用户登录，异步被调用

    Parameters
    ----------
    api_key : str
        用户的api_key
    """
    await asyncio.sleep(5)
    return


async def login(api_key: str = None):
    """用户登录

    Parameters
    ----------
    api_key : str
        用户的api_key
    """
    login_task = asyncio.create_task(_login(api_key))
    prefix = FONT.bold(FONT.blue("swanlab: "))
    loading_task = asyncio.create_task(FONT.loading("login...", interval=0.5, prefix=prefix))
    data = await login_task
    # 取消加载动画任务
    loading_task.cancel()
    # 最后需要刷去当前行，替换为新的内容，这将在外部完成
