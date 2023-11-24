#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 17:53:52
@File: test/start_server.py
@IDE: vscode
@Description:
        开启vue中编译后的后端服务，用于测试
"""

import time
from swanlab import SwanWeb

if __name__ == "__main__":
    sw = SwanWeb()
    print("start server")
    sw.run()
    # while True:
    #     print("running")
    #     time.sleep(2)
