#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/11 22:06
@File: pytest_watch.py
@IDE: pycharm
@Description:
    测试cli的watch命令
"""
from swanlab.cli.commands.dashboard.watch import get_free_port
import socket


def test_get_free_port():
    """
    测试get_free_port函数是否能正确获取一个可用端口
    """
    port = get_free_port()
    assert isinstance(port, int)
    assert port == 5092
    # 占用端口
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('0.0.0.0', 5092))
        port = get_free_port()
        assert isinstance(port, int)
        assert port != 5092
