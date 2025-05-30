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
import getpass
import sys

import requests
from swankit.env import is_windows
from swankit.log import FONT

from swanlab.api.info import LoginInfo
from swanlab.env import in_jupyter
from swanlab.error import ValidationError, APIKeyFormatError
from swanlab.log import swanlog
from swanlab.package import get_setting_url, get_host_api, get_host_web


def login_request(api_key: str, api_host: str, timeout: int = 20) -> requests.Response:
    """用户登录，请求后端接口完成验证"""
    resp = requests.post(url=f"{api_host}/login/api_key", headers={"authorization": api_key}, timeout=timeout)
    return resp


def login_by_key(api_key: str, timeout: int = 20, save: bool = True) -> LoginInfo:
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
    api_host, web_host = get_host_api(), get_host_web()
    resp = login_request(api_key, api_host, timeout)
    # api key写入token文件
    login_info = LoginInfo(resp, api_key, api_host, web_host)
    save and not login_info.is_fail and login_info.save()
    return login_info


def input_api_key(
    tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL + C' to quit",
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
        swanlog.info("You can find your API key at: " + FONT.yellow(get_setting_url()))
    # windows 额外打印提示信息
    if is_windows():
        tip += (
            '\nOn Windows, '
            f'use {FONT.yellow("Ctrl + Shift + V")} or {FONT.yellow("right-click")} '
            f'to paste the API key'
        )
    tip += ': '
    tip = FONT.swanlab(tip)
    ij = in_jupyter()
    ij and print(tip)
    key = getpass.getpass("" if ij else tip)
    sys.excepthook = _t
    return key


def code_login(api_key: str, save_key: bool = True) -> LoginInfo:
    """
    代码内登录，此时会覆盖本地token文件（非task模式下）
    :param api_key: 用户的api_key
    :param save_key: 是否保存api_key到本地token文件
    :return: 登录信息
    :raise ValidationError: 登录失败
    :raise APIKeyFormatError: api_key格式错误
    """
    APIKeyFormatError.check(api_key)
    tip = "Waiting for the swanlab cloud response."
    login_info: LoginInfo = FONT.loading(tip, login_by_key, args=(api_key, 20, save_key), interval=0.5)
    if login_info.is_fail:
        raise ValidationError("Login failed: " + str(login_info))
    return login_info


def terminal_login(api_key: str = None, save_key: bool = True) -> LoginInfo:
    """
    终端登录，此时直接覆盖本地token文件，但是新增交互，让用户输入api_key
    运行此函数，如果是认证失败的错误，重新要求用户输入api_key
    本地文件上层文件夹不保证存在，需要上层函数保证
    """
    # 1. api_key存在，跳过输入环节，直接请求登录接口，这与代码内swanlab.login方法一致
    # 2. api_key为None，提示用户输入
    is_input_key = api_key is None
    api_key = input_api_key() if api_key is None else api_key
    # 3. api_key校验失败且是输入的api_key，重新输入api_key，否则raise
    # 4. 登录失败且是输入的api_key，提示重新输入api_key，否则raise

    def login_again(error: Exception):
        swanlog.error(error)
        return input_api_key("Please try again, or press 'CTRL-C' to quit", True)

    while True:
        try:
            return code_login(api_key, save_key)
        # 登录失败且是输入的api_key，重新输入api_key
        except (APIKeyFormatError, ValidationError) as e:
            if not is_input_key:
                raise e
            api_key = login_again(e)


def _abort_tip(tp, _, __):
    """处理用户在input_api_key输入时按下CTRL+C的情况"""
    if tp == KeyboardInterrupt:
        sys.exit(0)
    # 如果不是CTRL+C，交给默认的异常处理
