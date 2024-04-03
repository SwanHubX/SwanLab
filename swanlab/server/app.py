#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 20:47:18
@File: swanlab\server\route.py
@IDE: vscode
@Description:
    综合服务 api
"""

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from swanlab.package import get_package_version
from .middleware.common import (
    resp_base,
    resp_static,
    catch_error,
    log_print,
    resp_params,
)

# 响应路径
from .settings import ASSETS

# 服务全局对象
app = FastAPI()
version = get_package_version()

# 注册前端静态文件路径
static_path = "/assets"
static = StaticFiles(directory=ASSETS)
app.mount(static_path, static)

# 将uvicorn的日志输出handler删除
import logging

# 删除 uvicorn logger
uvicorn_error = logging.getLogger("uvicorn.error")
uvicorn_error.disabled = True
uvicorn_access = logging.getLogger("uvicorn.access")
uvicorn_access.disabled = True


# ---------------------------------- 在此处注册中间件 ----------------------------------


@app.middleware("http")
async def _(*args, **kwargs):
    """基础中间件，调整响应结果，添加处理时间等信息"""
    return await resp_base(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    """资源中间件，此时所有与api相关的内容不会在此中间件中处理"""
    return await resp_static(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    """异常中间件，捕获异常，重构异常信息"""
    return await catch_error(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    """日志打印中间件"""
    return await log_print(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    """参数中间件，处理api请求中的参数校验问题，重新结构化校验错误结果

    参数校验错误并不会影响其他情况的响应结果
    此外由于参数校验错误在绝大多数情况应该是开发时的错误
    所以不会影响正式版本的性能
    """
    return await resp_params(*args, **kwargs)


# ---------------------------------- 在此处注册相关路由 ----------------------------------

# 导入数据相关的路由
from .router.experiment import router as experiment
from .router.project import router as project
from .router.namespace import router as namespace
from .router.chart import router as chart

# 媒体文件路由，允许前端获取其他产生的媒体文件
from .router.media import router as media

# 使用配置列表，统一导入
prefix = "/api/v1"
app.include_router(project, prefix=prefix + "/project")
app.include_router(experiment, prefix=prefix + "/experiment")
app.include_router(media, prefix=prefix + "/media")
app.include_router(namespace, prefix=prefix + "/namespace")
app.include_router(chart, prefix=prefix + "/chart")
