#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-19 19:00:37
@File: swanlab/utils/package.py
@IDE: vscode
@Description:
    用于管理swanlab的包管理器的模块，做一些封装
"""
from .env import get_package_path, get_save_dir
from .error import KeyFileError
from typing import Optional
import requests
import netrc
import json
import os


# ---------------------------------- 版本号相关 ----------------------------------

def get_package_version() -> str:
    """获取swanlab的版本号
    :return: swanlab的版本号
    """
    package_path = get_package_path()
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
    package_path = get_package_path()
    with open(package_path, "r", encoding="utf-8") as f:
        return json.load(f)["host"]["web"]


def get_host_api() -> str:
    """获取swanlab网站api网址
    :return: swanlab网站的api网址
    """
    package_path = get_package_path()
    with open(package_path, "r", encoding="utf-8") as f:
        return json.load(f)["host"]["api"]


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
    host = get_host_web()
    if not os.path.exists(path):
        raise KeyFileError("The file does not exist")
    nrc = netrc.netrc(path)
    info = nrc.authenticators(host)
    if info is None:
        raise KeyFileError(f"The host {host} does not exist")
    return info[2]


def save_key(username: str, password: str):
    """
    保存key到对应的文件目录下，文件名称为.netrc（basename）
    :param username: 保存的用户名
    :param password: 保存的密码
    :raises KeyFileError 传入的path路径文件名称不是.netrc或上级文件夹不存在
    """
    path = os.path.join(get_save_dir(), ".netrc")
    host = get_host_web()
    # 如果文件不存在，自动创建
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("")
    nrc = netrc.netrc(path)
    nrc.hosts[host] = (username, None, password)
    with open(path, "w") as f:
        f.write(nrc.__repr__())


def is_login() -> bool:
    """判断是否已经登录
    :return: 是否已经登录
    """
    try:
        _ = get_key()
        return True
    except KeyFileError:
        return False
