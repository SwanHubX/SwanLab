#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 21:12:52
@File: test/start_server.py
@IDE: vscode
@Description:
    开启后端测试服务，访问端口
    事实上在使用中并不是这样的，而是在命令行执行命令完成的，这里只是为了测试
"""

from swanlab import swanweb as sw
import sys
import os

print(os.getcwd())

print(sys.path)
import time

sw.init()


if __name__ == "__main__":
    print("start server")
    while True:
        print("running")
        time.sleep(2)
