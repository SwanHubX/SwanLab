"""
@author: cunyue
@file: session.py
@time: 2025/6/19 14:52
@description: 使用api_key登录后获取会话凭证
"""

import getpass
import sys
from typing import Union

import requests
from rich.status import Status
from rich.text import Text

from swanlab.env import is_windows, is_interactive
from swanlab.error import ValidationError, APIKeyFormatError, KeyFileError
from swanlab.log import swanlog
from swanlab.package import get_setting_url, get_host_api, get_host_web, fmt_web_host, save_key as sk, get_key
from ...session import create_session


class LoginInfo:
    """
    登录信息类，负责解析登录接口返回的信息，并且进行保存
    无论接口请求成功还是失败，都会初始化一个LoginInfo对象
    """

    def __init__(self, resp: requests.Response, api_key: str, api_host: str, web_host: str):
        self.__api_key = api_key
        self.__resp = resp
        self.__body = resp.json() if resp.status_code == 200 else {}
        self.__username = None
        """
        如果此属性不为None，username返回此属性
        """
        self.web_host = fmt_web_host(web_host)
        self.api_host = api_host

    @property
    def sid(self) -> Union[str, None]:
        """
        获取sid，如果请求失败则返回None
        """
        return self.__body.get("sid")

    @property
    def expired_at(self) -> Union[str, None]:
        """
        获取过期时间，如果请求失败则返回None
        """
        return self.__body.get("expiredAt")

    @property
    def username(self) -> Union[str, None]:
        """
        获取用户名，如果请求失败则返回None
        """
        if self.__username is not None:
            return self.__username
        return self.__body.get("userInfo", {}).get("username")

    @username.setter
    def username(self, value):
        """
        设置用户名
        """
        self.__username = value

    @property
    def is_fail(self):
        """
        判断登录是否失败
        """
        return self.__resp.status_code != 200

    @property
    def api_key(self):
        """
        获取api_key
        """
        if self.is_fail:
            return None
        return self.__api_key

    def __str__(self) -> str:
        """错误时会返回错误信息"""
        if self.__resp.reason == "OK":
            return "Login success"
        if self.__resp.reason == "Unauthorized" or self.__resp.reason == "Authorization Required":
            return "Error api key"
        if self.__resp.reason == "Forbidden":
            return "You need to be verified first"
        return str(self.__resp.status_code) + " " + self.__resp.reason

    def save(self):
        """
        保存登录信息
        """
        return sk(self.web_host, self.api_key, self.api_host)


def login_request(api_key: str, api_host: str, timeout: int = 20) -> requests.Response:
    """用户登录，请求后端接口完成验证"""
    session = create_session()
    resp = session.post(url=f"{api_host}/login/api_key", headers={"authorization": api_key}, timeout=timeout)
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
        swanlog.info("You can find your API key at:", Text(get_setting_url(), "yellow"))
    swanlog.info(tip, end='')
    # windows 额外打印提示信息
    if is_windows():
        swanlog.console.print(
            '\nOn Windows, ',
            'use',
            Text("Ctrl + Shift + V", 'yellow'),
            'or',
            Text("right-click", 'yellow'),
            'to paste the API key',
            end='',
        )
    swanlog.console.print(': ', end='')
    key = getpass.getpass("")
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
    with Status("Waiting for the swanlab cloud response.", spinner="dots"):
        login_info = login_by_key(api_key, 20, save_key)
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


def create_login_info(save: bool = True):
    """
    在代码运行时发起登录，获取登录信息，执行此方法会覆盖原有的login_info
    """
    key = None
    try:
        key = get_key()
    except KeyFileError:
        pass
    if key is None and not is_interactive():
        raise KeyFileError(
            "api key not configured (no-tty), call `swanlab.login(api_key=[your_api_key])` or set `swanlab.init(mode=\"local\")`."
        )
    return terminal_login(key, save)


def _abort_tip(tp, _, __):
    """处理用户在input_api_key输入时按下CTRL+C的情况"""
    if tp == KeyboardInterrupt:
        sys.exit(0)
    # 如果不是CTRL+C，交给默认的异常处理


__all__ = [
    "login_request",
    "login_by_key",
    "terminal_login",
    "code_login",
    "LoginInfo",
    "create_login_info",
]
