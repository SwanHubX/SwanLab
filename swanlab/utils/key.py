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
    """
    nrc = netrc.netrc(path)
    try:
        return nrc.authenticators(host)
    except Exception as e:
        raise KeyFileError("Failed to read token file") from e


def save_key(path: str, host: str, username: str, password: str):
    """
    保存key到对应的文件
    :param path: 保存位置
    :param host: 保存的host
    :param username: 保存的用户名
    :param password: 保存的密码
    :return:
    """
    # 如果文件不存在，自动创建
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("")
    nrc = netrc.netrc(path)
    nrc.hosts[host] = (username, None, password)
    with open(path, "w") as f:
        f.write(nrc.__repr__())
