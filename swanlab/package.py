#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-19 19:00:37
@File: swanlab/utils/package.py
@IDE: vscode
@Description:
    用于管理swanlab的包管理器的模块，做一些封装
"""
import json
import netrc
import os
import re
from typing import Optional

import requests

from .env import get_save_dir, SwanLabEnv
from .error import KeyFileError

package_path = os.path.join(os.path.dirname(__file__), "package.json")


# ---------------------------------- 版本号相关 ----------------------------------


def get_package_version() -> str:
    """获取swanlab的版本号
    :return: swanlab的版本号
    """
    # 读取package.json文件
    with open(package_path, "r") as f:
        return json.load(f)["version"]


def get_package_latest_version(timeout=0.5) -> Optional[str]:
    """
    获取swanlab的最新版本号
    :param timeout: 请求超时时间
    :return: 最新版本号
    """
    url = "https://pypi.org/pypi/swanlab/json"
    # noinspection PyBroadException
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            data = response.json()
            return data["info"]["version"]
        else:
            return None
    except Exception:
        return None


# ---------------------------------- 云端相关 ----------------------------------


class HostFormatter:
    def __init__(self, host: str = None, web_host: str = None):
        # 更新后的正则模式，允许匹配标准域名、IP地址和localhost
        self.pattern = re.compile(
            r'^(?:(https?)://)?'  # 可选协议 http 或 https
            r'('  # 主机部分（域名、IP、localhost）
            r'([a-zA-Z0-9.-]+\.[a-zA-Z]{2,63})'  # 标准域名
            r'|'  # 或
            r'(\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # IPv4地址
            r'|'  # 或
            r'(localhost)'  # 本地主机
            r')'
            r'(?::(\d{1,5}))?$'  # 可选端口号（1~5位数字）
        )
        self.host = host
        self.web_host = web_host

    def fmt(self, input_str: str) -> str:
        match = self.pattern.match(input_str.rstrip("/"))
        if match:
            protocol = match.group(1) or "https"  # 默认协议为 https
            host = match.group(2)
            port = match.group(6)

            # 构建标准化的 URL 输出
            result = f"{protocol}://{host}"
            if port:
                result += f":{port}"
            return result
        else:
            raise ValueError("Invalid host format")

    def __call__(self):
        """
        如果host或web_host不为空，格式化并设置环境变量
        :raises ValueError: host或web_host格式不正确
        """
        if self.host:
            try:
                os.environ[SwanLabEnv.API_HOST.value] = self.fmt(self.host) + "/api"
            except ValueError:
                raise ValueError("Invalid host: {}".format(self.host))
            self.web_host = self.host if self.web_host is None else self.web_host
        if self.web_host:
            try:
                os.environ[SwanLabEnv.WEB_HOST.value] = self.fmt(self.web_host)
            except ValueError:
                raise ValueError("Invalid web_host: {}".format(self.web_host))


def get_host_web() -> str:
    """获取swanlab网站网址
    :return: swanlab网站的网址
    """
    return os.getenv(SwanLabEnv.WEB_HOST.value, "https://swanlab.cn")


def get_host_api() -> str:
    """获取swanlab网站api网址
    :return: swanlab网站的api网址
    """
    return os.getenv(SwanLabEnv.API_HOST.value, "https://api.swanlab.cn/api")


def fmt_web_host(web_host: str = None) -> str:
    """
    如果web_host为None，则使用默认的web_host
    并且格式化web_host，去除结尾的/
    :param web_host: web_host
    :return: 格式化后的web_host
    """
    if web_host is None:
        web_host = get_host_web()
    return web_host.rstrip("/")


def get_setting_url(web_host: str = None) -> str:
    """获取用户设置的url
    与实验相关的url不同，这个url在http对象之前被使用，因此不绑定在http对象中
    :return: 用户设置的url
    """
    return fmt_web_host(web_host) + "/space/~/settings"


def get_login_url(web_host: str = None) -> str:
    """获取登录的url
    与实验相关的url不同，这个url在http对象之前被使用，因此不绑定在http对象中
    :return: 登录的url
    """
    return fmt_web_host(web_host) + "/login"


# ---------------------------------- 登录相关 ----------------------------------


def get_nrc_path() -> str:
    """
    获取netrc文件路径
    """
    return os.path.join(get_save_dir(), ".netrc")


def get_key():
    """使用标准netrc库解析token文件，获取token
    :raise KeyFileError: 文件不存在或者host不存在
    :return: token
    """
    env_key = os.getenv(SwanLabEnv.API_KEY.value)
    if env_key is not None:
        return env_key
    path, host = get_nrc_path(), get_host_api().rstrip("/api")
    if not os.path.exists(path):
        raise KeyFileError("The file does not exist")
    nrc = netrc.netrc(path)
    info = nrc.authenticators(host)

    # 向下兼容 https://github.com/SwanHubX/SwanLab/issues/792#issuecomment-2649685881
    if info is None:
        info = nrc.authenticators(host + "/api")
        if info is not None:
            nrc.hosts = {host: info}
            with open(path, "w") as f:
                f.write(repr(nrc))

    if info is None:
        raise KeyFileError(f"The host {host} does not exist")
    return info[2]


def save_key(username: str, password: str, host: str = None) -> bool:
    """
    保存key到对应的文件目录下，文件名称为.netrc（basename）
    此函数不考虑上层文件存在的清空，但是会在调用的get_save_dir()函数中进行检查
    :param username: 保存的用户名，默认为user，可选择存储为前端网页ip或者域名
    :param password: 保存的密码
    :param host: 保存的host
    :return: 是否保存，如果已经存在，则不保存
    """
    if host is None:
        host = get_host_api()
    host = host.rstrip("/api")
    path = get_nrc_path()
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("")
    nrc = netrc.netrc(path)
    new_info = (username, "", password)
    # 避免重复的写
    info = nrc.authenticators(host)
    if info is None or (info[0], info[2]) != (new_info[0], new_info[2]):
        # 同时只允许存在一个host： https://github.com/SwanHubX/SwanLab/issues/797
        nrc.hosts = {host: new_info}
        with open(path, "w") as f:
            f.write(nrc.__repr__())
        return True
    return False


class LoginCheckContext:
    """
    进入上下文时，会删除环境变量中的api key，退出上下文时会恢复原来的值
    """

    def __init__(self):
        self.__tmp_key = None
        """
        临时保存的key
        """
        self.is_login = False
        """
        标注是否已经登录
        """

    def __enter__(self):
        self.__tmp_key = os.environ.get(SwanLabEnv.API_KEY.value)
        if self.__tmp_key is not None:
            del os.environ[SwanLabEnv.API_KEY.value]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # 恢复原来的值
        if self.__tmp_key is not None:
            os.environ[SwanLabEnv.API_KEY.value] = self.__tmp_key
        if exc_type is KeyFileError:  # 未登录
            return True
        elif exc_type is not None:  # 其他错误
            return False
        self.is_login = True
        return True


def has_api_key() -> bool:
    """判断是否已经登录，与当前的host相关
    如果环境变量中有api key，则认为已经登录
    但不会检查key的有效性
    :return: 是否已经登录
    """
    with LoginCheckContext() as checker:
        _ = get_key()
    return checker.is_login
