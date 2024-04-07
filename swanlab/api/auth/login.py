#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 15:30:24
@File: swanlab/auth/login.py
@IDE: vscode
@Description:
    用户登录接口，输入用户的apikey，保存用户token到本地
    进行一些交互定义和数据请求
"""
import asyncio
from swanlab.error import ValidationError
from swanlab.utils import FONT
from swanlab.package import get_user_setting_path, get_host_api
from .info import LoginInfo
import getpass
import sys
import requests
import time


def login_request(api_key: str, timeout: int = 20) -> requests.Response:
    """用户登录，请求后端接口完成验证"""
    return requests.post(f"{get_host_api()}/login/api_key", headers={'authorization': api_key}, timeout=timeout)


async def login_by_key(api_key: str, timeout: int = 20, save: bool = True) -> LoginInfo:
    """用户登录，异步调用接口完成验证
    返回后端内容(dict)，如果后端请求失败，返回None

    Parameters
    ----------
    api_key : str
        用户api_key
    timeout : int, optional
        请求认证的超时时间，单位秒
    save : bool, optional
        是否保存到本地token文件
    """
    try:
        resp = login_request(api_key, timeout)
    except requests.exceptions.RequestException:
        # 请求超时等网络错误
        raise ValidationError("Network error, please try again.")
    # api key写入token文件
    login_info = LoginInfo(resp, api_key)
    save and not login_info.is_fail and login_info.save()
    await asyncio.sleep(2)

    return login_info


def input_api_key(
    tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL-C' to quit: ",
    again: bool = False,
) -> str:
    """让用户输入apikey
    此时有两条消息，第一条消息为固定格式，第二条消息

    Parameters
    ----------
    tip : str
        提示信息
    again : bool, optional
        是否是重新输入api_key，如果是，不显示额外的提示信息
    """
    _t = sys.excepthook
    sys.excepthook = _abort_tip
    if not again:
        print(FONT.swanlab("Logging into swanlab cloud."))
        print(FONT.swanlab("You can find your API key at: " + get_user_setting_path()))
    key = getpass.getpass(FONT.swanlab(tip))
    sys.excepthook = _t
    return key


async def code_login(api_key: str) -> LoginInfo:
    """
    代码内登录，此时会覆盖本地token文件

    Parameters
    ----------
    api_key : str
        用户api_key
    """
    login_task = asyncio.create_task(login_by_key(api_key))
    tip = "Waiting for the swanlab cloud response."
    loading_task = asyncio.create_task(FONT.loading(tip, interval=0.5))
    login_info: LoginInfo = await login_task
    # 取消加载动画任务
    loading_task.cancel()
    # 最后需要刷去当前行, 不再显示加载动画
    FONT.brush("", length=100)
    if login_info.is_fail:
        print(FONT.swanlab("Login failed: " + str(login_info).lower(), color="red"))
        raise ValidationError("Login failed: " + str(login_info))
    return login_info


def terminal_login(api_key: str = None) -> LoginInfo:
    """
    终端登录，此时直接覆盖本地token文件，但是新增交互，让用户输入api_key
    运行此函数，如果是认证失败的错误，重新要求用户输入api_key
    """
    # 1. api_key存在，跳过输入环节，直接请求登录接口，这与代码内swanlab.login方法一致
    # 2. api_key为None，提示用户输入
    input_key = api_key is None
    if api_key is None:
        api_key = input_api_key()
    while True:
        try:
            return asyncio.run(code_login(api_key))
        # 如果是登录失败且是输入的api_key，提示重新输入api_key
        except ValidationError as e:
            if input_key:
                api_key = input_api_key("Please try again, or press 'CTRL-C' to quit: ", True)
            else:
                raise e


def _abort_tip(tp, val, tb):
    """处理用户在input_api_key输入时按下CTRL+C的情况"""
    if tp == KeyboardInterrupt:
        print("\n" + FONT.swanlab("Aborted!", color="red"))
        sys.exit(0)
    # 如果不是CTRL+C，交给默认的异常处理
