#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 15:30:24
@File: swanlab/auth/login.py
@IDE: vscode
@Description:
    用户登录接口，输入用户的apikey，保存用户token到本地
"""
import asyncio
from ..error import ValidationError
from ..utils import FONT
from ..utils.package import user_setting_path


async def _login(user_token: str):
    """用户登录，异步调用接口完成验证

    Parameters
    ----------
    user_token : str
        用户api_key
    """
    await asyncio.sleep(5)
    return


async def _check_key_format(api_key: str):
    """检查api_key的格式是否正确

    Parameters
    ----------
    api_key : str
        用户api_key
    """
    if len(api_key) < 10:
        raise ValidationError("api_key格式错误")
    return


def input_api_key(tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL-C' to quit: "):
    """让用户输入apikey
    此时有两条消息，第一条消息为固定格式，第二条消息

    Parameters
    ----------
    str : str
        用户api_key

    """
    print(FONT.swanlab("You can find your API key at: " + user_setting_path()))

    return "token"


def code_login():
    """
    代码内登录，此时会覆盖本地token文件
    """


def terminal_login(api_key: str = None):
    """
    终端登录，此时直接覆盖本地token文件，但是新增交互，让用户输入api_key
    """
    # 1. api_key存在，跳过输入环节，直接请求登录接口，这与代码内swanlab.login方法一致
    # 2. api_key为None，提示用户输入
    if api_key is None:
        api_key = input_api_key("Please input your api_key: ")
