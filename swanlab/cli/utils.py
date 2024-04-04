#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-15 15:53:21
@File: swanlab/cli/utils.py
@IDE: vscode
@Description:
    命令行工具
"""
import os
import psutil
import socket
import click
from ..utils import FONT, file
from ..env import PORT, HOST, ROOT


def is_valid_ip(ctx, param, ip: str) -> None:
    """检测输入的是否是合法的ip地址,完成环境变量的注入

    Parameters
    ----------
    ctx : click.Context
        上下文
    param : click.Parameter
        参数
    ip : str
        带检测的字符串
    """
    if ip is None:
        return
    if not file.is_ipv4(ip):
        raise click.BadParameter("Invalid ip address: " + ip)
    os.environ[HOST] = ip


def is_valid_port(ctx, param, port: int) -> int:
    """检测是否是合法的端口号

    Parameters
    ----------
    ctx : click.Context
        上下文
    param : click.Parameter
        参数
    port : int
        带检测的端口号
    """
    if port is None:
        return
    if not file.is_port(port):
        raise click.BadParameter("Invalid port number: " + str(port))
    os.environ[PORT] = str(port)


def is_valid_root_dir(ctx, param, log_dir: str) -> str:
    """检测是否是合法的日志目录，保证其可读且存在

    Parameters
    ----------
    ctx : click.Context
        上下文
    param : click.Parameter
        参数
    log_dir : str
        带检测的日志目录
    """
    # 将日志目录注入环境变量，在这之前先转换为绝对路径

    if log_dir is None:
        return

    # 将传入的路径转换为绝对路径
    log_dir = os.path.abspath(log_dir)

    # 必须是一个绝对路径
    if not os.path.isabs(log_dir):
        raise click.BadParameter("Log dir must be an absolute path: " + log_dir)
    # 路径必须存在
    if not os.path.isdir(log_dir):
        raise click.BadParameter("Log dir is not a directory: " + log_dir)
    # 路径必须可读
    if not os.access(log_dir, os.R_OK):
        raise click.BadParameter("Log dir is not readable: " + log_dir)

    os.environ[ROOT] = log_dir


class URL(object):
    # 生成链接提示,先生成各个组件
    _arrow = "\t\t\t➜"
    arrow = FONT.bold(FONT.green(_arrow))
    local = arrow + FONT.bold("  Local:   ")
    netwo = arrow + FONT.bold("  Network: ")

    def __init__(self, ip, port) -> None:
        self.ip = ip
        self.port = port

    def __str__(self) -> str:
        url = FONT.blue(f"http://{self.ip}:{self.port}")
        if self.is_localhost(self.ip):
            return self.local + url
        else:
            return self.netwo + url

    @classmethod
    def last_tip(cls) -> str:
        """
        打印最后一条提示信息
        """
        t = FONT.dark_gray("  press ") + FONT.bold(FONT.default("ctrl + c")) + FONT.dark_gray(" to quit")
        return FONT.dark_green(cls._arrow) + t

    @staticmethod
    def is_localhost(ip):
        return ip == "127.0.0.1" or ip == "localhost"

    @staticmethod
    def is_zero_ip(ip):
        return ip == "0.0.0.0"

    @staticmethod
    def get_all_ip() -> list:
        """获取所有可用的ip地址

        Parameters
        ----------
        ip : str
            本机的ip地址

        Returns
        -------
        tuple
            所有可用的ip地址
        """
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
        # 对 IPv4 进行排序，"127.0.0.1" 在最前面，剩下按照从小到大排序
        ipv4.sort(key=lambda x: (x != "127.0.0.1", x), reverse=True)
        ipv4.reverse()
        return ipv4
