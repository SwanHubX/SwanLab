#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:35:12
@File: swanlab/utils/key.py
@IDE: vscode
@Description:
    token文件操作，不作为模块方法暴露，如果需要使用相关方法需要指定此文件
"""
from ..error import KeyFileError
import netrc
import os


def get_key(path: str, host: str):
    """使用标准netrc库解析token文件，获取token

    Parameters
    ----------
    path : str
        token文件路径
    host : str
        token对应的host

    Raises
    ------
    KeyFileError
        传入的path路径文件不存在
    """
    if not os.path.exists(path):
        raise KeyFileError("The file does not exist")
    nrc = netrc.netrc(path)
    info = nrc.authenticators(host)
    if info is None:
        raise KeyFileError(f"The host {host} does not exist")
    return info


def save_key(path: str, host: str, username: str, password: str):
    """
    保存key到对应的文件目录下，文件名称为.netrc（basename）
    :param path: 保存位置，必须是文件夹
    :param host: 保存的host
    :param username: 保存的用户名
    :param password: 保存的密码
    :raises KeyFileError 传入的path路径文件名称不是.netrc或上级文件夹不存在
    """
    # 传入的path路径文件名称不是.netrc
    if os.path.basename(path) != ".netrc":
        raise KeyFileError("The file name must be .netrc")
    # 上级文件夹不存在
    if not os.path.exists(os.path.dirname(path)):
        raise KeyFileError("The parent folder does not exist")
    # 如果文件不存在，自动创建
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("")
    nrc = netrc.netrc(path)
    nrc.hosts[host] = (username, None, password)
    with open(path, "w") as f:
        f.write(nrc.__repr__())
