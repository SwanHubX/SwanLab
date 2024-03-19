#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:22:12
@File: swanlab/auth/experiment.py
@IDE: vscode
@Description:
    实验认证接口，此时用户已登录，但是可能认证失败，需要重新认证
"""
import asyncio
import requests
from ..utils import FONT


async def _get_exp_token(user_token: str):
    """用户登录，异步调用

    Parameters
    ----------
    user_token : str
        用户api_key经过后端验证后的token
    """
    await asyncio.sleep(5)
    return


async def get_exp_token(user_token: str = None):
    """通过apikey获取实验令牌
    接下来通过此令牌上传实验日志
    获取令牌的途中显示转圈圈，表示正在获取

    Parameters
    ----------
    user_token : str
        用户的token
    """
    login_task = asyncio.create_task(_get_exp_token(user_token))
    prefix = FONT.bold(FONT.blue("swanlab: "))
    loading_task = asyncio.create_task(FONT.loading("login...", interval=0.5, prefix=prefix))
    data = await login_task
    # 取消加载动画任务
    loading_task.cancel()
    # 最后需要刷去当前行, 不再显示加载动画
    FONT.brush("")
    return "token"
