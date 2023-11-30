#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 21:37:26
@File: swanlab\server\router.py
@IDE: vscode
@Description:
    服务路由模块，用于实现服务的路由功能
"""

from fastapi import FastAPI
from fastapi.responses import HTMLResponse, FileResponse

# 响应路径
from ..env import INDEX, LOGO

# 导入数据相关的路由
from .api.data import router as data_router


# 服务全局对象
app = FastAPI()


# 响应首页
@app.get("/", response_class=HTMLResponse)
async def _():
    # 读取 HTML 文件内容并返回
    with open(INDEX, "r", encoding="utf-8") as file:
        html_content = file.read()
    return HTMLResponse(content=html_content, status_code=200)


# 响应logo内容
# TODO 后续可以考虑将logo.ico放在assets中，这样就不需要单独响应了
@app.get("/logo.ico")
async def _():
    return FileResponse(LOGO)


app.include_router(data_router)
