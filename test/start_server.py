#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-24 21:12:52
@File: test/start_server.py
@IDE: vscode
@Description:
    开启后端测试服务，访问端口
    事实上在使用中并不是这样的，而是在命令行执行命令完成的，这里只是为了测试和开发
    增加了实际开发中不会用到的热启动功能
"""
from swanlab.server.router import app
import uvicorn


if __name__ == "__main__":
    uvicorn.run("start_server:app", host="0.0.0.0", port=6092, reload=True, log_level="critical")
    # swl.info("hello")
