#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 21:12:52
@File: test/start_server.py
@IDE: vscode
@Description:
    开启后端测试服务，访问端口
"""
from swanlab import SwanWeb
import time

if __name__ == "__main__":
    sw = SwanWeb()
    print("start server")
    sw.run()
    while True:
        print("running")
        time.sleep(2)
