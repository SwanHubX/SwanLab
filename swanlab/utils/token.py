#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 19:35:12
@File: swanlab/utils/token.py
@IDE: vscode
@Description:
    token文件操作，不作为模块方法暴露，如果需要使用相关方法需要指定此文件
"""
from ..error import TokenFileError
import netrc


def get_token(path: str, host: str):
    """使用标准netrc库解析token文件，获取token

    Parameters
    ----------
    path : str
        token文件路径
    """
    nrc = netrc.netrc(path)
    try:
        return nrc.authenticators(host)
    except Exception as e:
        raise TokenFileError("Failed to read token file") from e


def save_token(path: str, host: str, username: str, password: str):
    """保存token到token文件"""
    nrc = netrc.netrc()
    nrc.hosts[host] = (username, None, password)
    with open(path, "w") as f:
        f.write(nrc.__repr__())
