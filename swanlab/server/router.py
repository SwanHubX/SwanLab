#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 20:47:18
@File: swanlab\server\route.py
@IDE: vscode
@Description:
    综合服务 api
"""

from fastapi import FastAPI, status
from fastapi.responses import HTMLResponse, FileResponse, Response
from fastapi.staticfiles import StaticFiles
import ujson
import time

# 响应路径
from ..env import INDEX, LOGO, ASSETS


# 服务全局对象
app = FastAPI()

# 注册静态文件路径
static_path = "/assets"
static = StaticFiles(directory=ASSETS)
app.mount(static_path, static)


# ---------------------------------- 在此处注册中间件 ----------------------------------


@app.middleware("http")
async def _(request, call_next):
    """基础中间件，调整响应结果，添加处理时间等信息"""
    # 如果请求路径不以'/api'开头，说明并不是后端服务的请求，直接返回
    if not request.url.path.startswith("/api"):
        return await call_next(request)
    # 记录请求开始时间
    start_time = time.time()
    # 调用下一个中间件或者最终的路由处理函数
    # 我们约定路由处理函数最终返回三个参数，一个是错误码，一个是错误信息，一个是响应内容
    response = await call_next(request)
    # 记录请求结束时间
    end_time = time.time()
    # 计算处理时间，添加到响应头中
    process_time = round(end_time - start_time, 4)
    response.headers["X-Process-Time"] = str(process_time)
    # 返回响应
    return response


# ---------------------------------- 在此处注册相关路由 ----------------------------------


# 导入数据相关的路由
from .api.test import router as test
from .api.project import router as project
from .api.experiment import router as experiment


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


# ---------------------------------- 加载动态路由 ----------------------------------


# TODO 使用配置列表，统一导入
prefix = "/api/v1"
app.include_router(test, prefix=prefix)
app.include_router(project, prefix=prefix + "/project")
app.include_router(experiment, prefix=prefix + "/experiments")
