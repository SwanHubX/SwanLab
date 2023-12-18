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
import click


def is_vaild_ip(ctx, param, ip: str) -> tuple:
    """检测是否是合法的ip地址

    Parameters
    ----------
    ctx : click.Context
        上下文
    param : click.Parameter
        参数
    ip : str
        带检测的字符串
    """
    ip = str(ip)
    pattern = re.compile(r"^((2[0-4]\d|25[0-5]|[01]?\d\d?)\.){3}(2[0-4]\d|25[0-5]|[01]?\d\d?)$")
    if not pattern.match(ip):
        raise click.BadParameter("Invalid IP address format: " + ip)
    # 没有问题，获取当前机器的所有ip地址
    interfaces = psutil.net_if_addrs()
    ipv4 = []
    # APIPA 地址范围
    apipa_range = range(169, 255)
    for _, addresses in interfaces.items():
        for address in addresses:
            # 如果是ipv4地址，且可以被访问到
            if address.family == socket.AddressFamily.AF_INET:
                # 排除 APIPA 地址范围
                octets = list(map(int, address.address.split(".")))
                if octets[0] == 169 and octets[1] in apipa_range:
                    continue
                ipv4.append(address.address)
    if ip not in ipv4 and ip != "0.0.0.0":
        raise click.BadParameter("IP address '" + ip + "' should be one of " + str(ipv4) + ".")
    return ip, ipv4


def is_available_port(host, port):
    """检测端口是否可用

    Parameters
    ----------
    host : str
        ip地址
    port : int
        端口号
    """
    try:
        with socket.create_server((host, port)):
            return True
    except:
        return False
