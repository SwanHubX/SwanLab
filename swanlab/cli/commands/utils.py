#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:27
@File: utils.py
@IDE: pycharm
@Description:
    
"""
import sys
from typing import Optional
from swanlab.api import get_http
from swanlab.error import ApiError
from swanlab.log import swanlog


class UseTaskHttp:
    """
    主要用于检测http响应是否为3xx字段，如果是则要求用户更新版本
    使用此类之前需要先调用login_init_sid()函数完成全局http对象的初始化
    """

    def __init__(self):
        self.http = get_http()

    def __enter__(self):
        return self.http

    def __exit__(self, exc_type, exc_val: Optional[ApiError], exc_tb):
        if exc_type is ApiError:
            # api已过期，需要更新swanlab版本
            if exc_val.resp.status_code // 100 == 3:
                swanlog.info("SwanLab in your environment is outdated. Upgrade: `pip install -U swanlab`")
                sys.exit(3)
        return False
