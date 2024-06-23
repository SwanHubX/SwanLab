#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-19 19:00:37
@File: swanlab/utils/package.py
@IDE: vscode
@Description:
    用于管理swanlab的包管理器的模块，做一些封装
"""
from .env import get_save_dir, SwanLabEnv
from .error import KeyFileError
from typing import Optional
import requests
import netrc
import json
import os

package_path = os.path.join(os.path.dirname(__file__), "package.json")


# ---------------------------------- 版本号相关 ----------------------------------

def get_package_version() -> str:
    """获取swanlab的版本号
    :return: swanlab的版本号
    """
    if SwanLabEnv.SWANLAB_VERSION.value in os.environ:
        return os.environ[SwanLabEnv.SWANLAB_VERSION.value]
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


# ---------------------------------- 云端url相关 ----------------------------------


def get_host_web() -> str:
    """获取swanlab网站网址
    :return: swanlab网站的网址
    """
    return os.getenv(SwanLabEnv.SWANLAB_WEB_HOST.value, "https://swanlab.cn")


def get_host_api() -> str:
    """获取swanlab网站api网址
    :return: swanlab网站的api网址
    """
    return os.getenv(SwanLabEnv.SWANLAB_API_HOST.value, "https://swanlab.cn/api")


def get_user_setting_path() -> str:
    """获取用户设置的url
    :return: 用户设置的url
    """
    return get_host_web() + "/settings"


def get_project_url(username: str, projname: str) -> str:
    """获取项目的url
    :param username: 用户名
    :param projname: 项目名
    :return: 项目的url
    """
    return get_host_web() + "/" + username + "/" + projname


def get_experiment_url(username: str, projname: str, expid: str) -> str:
    """获取实验的url
    :param username: 用户名
    :param projname: 项目名
    :param expid: 实验id
    :return: 实验的url
    """
    return get_project_url(username, projname) + "/" + expid


# ---------------------------------- 登录相关 ----------------------------------


def get_key():
    """使用标准netrc库解析token文件，获取token
    :raise KeyFileError: 文件不存在或者host不存在
    :return: token
    """
    path = os.path.join(get_save_dir(), ".netrc")
    host = get_host_api()
    if not os.path.exists(path):
        raise KeyFileError("The file does not exist")
    nrc = netrc.netrc(path)
    info = nrc.authenticators(host)
    if info is None:
        raise KeyFileError(f"The host {host} does not exist")
    return info[2]


def save_key(username: str, password: str, host: str = None):
    """
    保存key到对应的文件目录下，文件名称为.netrc（basename）
    此函数不考虑上层文件存在的清空，但是会在调用的get_save_dir()函数中进行检查
    :param username: 保存的用户名
    :param password: 保存的密码
    :param host: 保存的host
    """
    if host is None:
        host = get_host_api()
    path = os.path.join(get_save_dir(), ".netrc")
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("")
    nrc = netrc.netrc(path)
    nrc.hosts[host] = (username, None, password)
    with open(path, "w") as f:
        f.write(nrc.__repr__())


def is_login() -> bool:
    """判断是否已经登录，与当前的host相关
    但不会检查key的有效性
    :return: 是否已经登录
    """
    try:
        _ = get_key()
        return True
    except KeyFileError:
        return False
