#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-19 22:10:15
@File: swanlab\server\middleware\common.py
@IDE: vscode
@Description:
    常规中间件
"""

from fastapi import Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
import time
from ..module.resp import UNEXPECTED_ERROR_500, PARAMS_ERROR_422
from ...log import swanlog as swl
from swanlab.package import get_package_version
from ..settings import ASSETS, INDEX  # 响应路径

version = get_package_version()

# 注册静态文件路径
static_path = "/assets"
static = StaticFiles(directory=ASSETS)


async def resp_base(request, call_next):
    """基础中间件，调整响应结果，添加处理时间等信息"""
    # 如果请求路径不以 '/api' 开头，说明并不是后端服务的请求，直接返回
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
    response.headers["SwanLab-Process-Time"] = str(process_time)
    # 添加版本信息
    response.headers["SwanLab-Version"] = version
    # 返回响应
    return response


async def resp_static(request, call_next):
    """资源中间件，此时所有与api相关的内容不会在此中间件中处理"""
    if request.url.path.startswith(static_path):
        # 如果是请求静态资源，直接返回
        return await call_next(request)
    if request.url.path.startswith("/api"):
        # 如果是请求api，直接返回
        return await call_next(request)
    # 剩余情况返回index.html，交由前端路由处理
    with open(INDEX, "r", encoding="utf-8") as file:
        html_content = file.read()
    return HTMLResponse(content=html_content, status_code=200)


async def catch_error(request: Request, call_next):
    """异常中间件，捕获异常，重构异常信息"""
    if not request.url.path.startswith("/api"):
        # 如果不是请求api，直接返回
        return await call_next(request)
    try:
        return await call_next(request)
    except Exception as e:
        swl.error(str(e))
        return UNEXPECTED_ERROR_500(str(e))


async def log_print(request: Request, call_next):
    """日志打印中间件"""
    swl.debug("[" + request.method + "] from " + request.base_url._url)
    resp = await call_next(request)
    # 拿到状态码
    # status = str(resp.status_code)
    if not request.url.path.startswith("/api"):
        # 如果不是请求api，直接返回
        swl.debug("[" + str(resp.status_code) + "] " + request.method + " assets: " + request.url.path)
    else:
        content = "[" + str(resp.status_code) + "] " + request.method + " api: " + request.url.path
        swl.debug(content)
    return resp


async def resp_params(request: Request, call_next):
    """参数中间件，处理api请求中的参数校验问题，重新结构化校验错误结果"""
    if not request.url.path.startswith("/api"):
        # 如果不是请求api，直接返回
        return await call_next(request)
    resp = await call_next(request)
    # 拿到状态码
    status = resp.status_code
    if status == 422:
        # 参数校验错误，重构错误信息
        # 拿到响应体
        import json

        body = [chunk async for chunk in resp.body_iterator][0].decode()
        body = json.loads(body)
        detail = body["detail"][0]
        msg = detail["msg"].split(" ")[1:]
        loc = detail["loc"][1]
        msg.insert(0, loc)
        msg = " ".join(msg)
        return PARAMS_ERROR_422(msg)
    """
    参数校验错误并不会影响其他情况的响应结果
    此外由于参数校验错误在绝大多数情况应该是开发时的错误
    所以不会影响正式版本的性能
    """
    return resp
