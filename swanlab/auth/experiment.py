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
from ..env import is_login, get_user_api_key
from ..error import NotLoginError, TokenFileError


async def _get_exp_token(user_token: str):
    """用户登录，异步调用

    Parameters
    ----------
    user_token : str
        用户api_key经过后端验证后的token
    """
    await asyncio.sleep(5)
    # TODO 如果token是12345，模拟错误
    if user_token == "12345":
        return None
    return "token"


async def get_exp_token():
    """通过apikey获取实验令牌
    接下来通过此令牌上传实验日志
    获取令牌的途中显示转圈圈，表示正在获取
    """
    if not is_login():
        raise NotLoginError("Please login first")
    # 此时get_user_api_key必然成功
    api_key = get_user_api_key()
    login_task = asyncio.create_task(_get_exp_token(api_key))
    # 显示加载动画
    prefix = FONT.bold(FONT.blue("swanlab: "))
    tip = "Creating experiment in swanlab cloud..."
    loading_task = asyncio.create_task(FONT.loading(tip, interval=0.5, prefix=prefix))
    data = await login_task
    # 取消加载动画任务
    loading_task.cancel()
    # 最后需要刷去当前行, 不再显示加载动画
    FONT.brush("")
    # 在此完成错误处理，比如后端请求失败之类的，直接抛出错误
    if data is None:
        raise TokenFileError("Failed to get experiment token: 500")
    # TODO 其他错误就直接返回状态码

    return "token"
