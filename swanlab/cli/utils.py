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
from ..utils import FONT


def is_vaild_ip(ctx, param, ip: str) -> tuple:
    """检测是否是合法的ip地址,最终会返回所有可访问的地址,并且排序
    127.0.0.1在最前面,剩下按照从小到大排序

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
    ipv4: list = []
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
    # 对 IPv4 进行排序，"127.0.0.1" 在最前面，剩下按照从小到大排序
    ipv4.sort(key=lambda x: (x != "127.0.0.1", x), reverse=True)
    ipv4.reverse()
    return ip, ipv4


class URL(object):
    # 生成链接提示,先生成各个组件
    arrow = FONT.green("\t\t\t➜")
    local = arrow + "  Local:   "
    netwo = arrow + "  Network: "

    def __init__(self, ip, port) -> None:
        self.ip = ip
        self.port = port

    def __str__(self) -> str:
        url = FONT.bold(f"http://{self.ip}:{self.port}")
        if self.is_localhost(self.ip):
            return self.local + url
        else:
            return self.netwo + url

    @staticmethod
    def is_localhost(ip):
        return ip == "127.0.0.1" or ip == "localhost"

    @staticmethod
    def is_zero_ip(ip):
        return ip == "0.0.0.0"
