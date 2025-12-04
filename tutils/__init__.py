#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/21 17:01
@File: __init__.py
@IDE: pycharm
@Description:
    tutils模块的初始化文件
"""
from .check import *
from .config import *


def reset_some_env():
    os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = SWANLOG_FOLDER
    os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = SWANLAB_FOLDER
    os.environ[SwanLabEnv.RUNTIME.value] = runtime
    SwanLabEnv.set_default()
    SwanLabEnv.check()
    if not is_skip_cloud_test:
        os.environ[SwanLabEnv.API_HOST.value] = API_HOST
        os.environ[SwanLabEnv.WEB_HOST.value] = WEB_HOST
        os.environ[SwanLabEnv.API_KEY.value] = API_KEY
    if SwanLabEnv.LOG_LEVEL.value in os.environ:
        del os.environ[SwanLabEnv.LOG_LEVEL.value]


# env 必须在 config 之后导入
from swanlab.env import SwanLabEnv

API_HOST = os.getenv(SwanLabEnv.API_HOST.value)
WEB_HOST = os.getenv(SwanLabEnv.WEB_HOST.value)
API_KEY = os.getenv(SwanLabEnv.API_KEY.value)


if not os.path.exists(TEMP_PATH):
    os.mkdir(TEMP_PATH)
reset_some_env()


def open_dev_mode() -> str:
    """
    开启开发模式，此时会返回开发环境的api-key并且创建测试目录
    在上层config部分已经执行了环境变量注入
    :return: api-key
    """
    return API_KEY
