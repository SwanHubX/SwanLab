#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 20:47:18
@File: swanlab\server\route.py
@IDE: vscode
@Description:
    综合服务 api
"""

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
import time
from .module.resp import UNEXPECTED_ERROR_500, PARAMS_ERROR_422
from ..log import swanlog as swl
from ..utils import get_package_version

from .middleware.common import (
    resp_base,
    resp_static,
    catch_error,
    log_print,
    resp_params,
)

# 响应路径
from .settings import ASSETS, INDEX


# 服务全局对象
app = FastAPI()
version = get_package_version()

# 注册静态文件路径
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
    return await resp_base(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    return await resp_static(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    return await catch_error(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    return await log_print(*args, **kwargs)


@app.middleware("http")
async def _(*args, **kwargs):
    return await resp_params(*args, **kwargs)


# ---------------------------------- 在此处注册相关路由 ----------------------------------


# 导入数据相关的路由
# from .api.project import router as project
# from .api.experiment import router as experiment
from .routes.experiment import router as experiment
from .routes.project import router as project

# 使用配置列表，统一导入
prefix = "/api/v1"
app.include_router(project, prefix=prefix + "/project")
app.include_router(experiment, prefix=prefix + "/experiment")
