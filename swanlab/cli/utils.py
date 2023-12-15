#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 15:53:21
@File: swanlab/cli/utils.py
@IDE: vscode
@Description:
    命令行工具
"""
import re
import psutil
import socket


def is_vaild_ip(ip: str):
    """检测是否是合法的ip地址

    Parameters
    ----------
    ip : str
        带检测的字符串
    """
    ip = str(ip)
    pattern = re.compile(r"^((2[0-4]\d|25[0-5]|[01]?\d\d?)\.){3}(2[0-4]\d|25[0-5]|[01]?\d\d?)$")
    if not pattern.match(ip):
        return False
    return True


def get_all_ip():
    """获取当前机器的所有可用ip地址"""
    interfaces = psutil.net_if_addrs()
    ipv4 = []
    for _, addresses in interfaces.items():
        for address in addresses:
            # 如果是ipv4地址
            if address.family == socket.AddressFamily.AF_INET:
                ipv4.append(address.address)
    return ipv4
